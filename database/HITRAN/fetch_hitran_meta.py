#!/usr/bin/env python3
"""Fetch HITRAN isotopologue metadata and partition function files."""

import os
import re
import requests
from concurrent.futures import ThreadPoolExecutor, as_completed
from tqdm import tqdm

URL = 'https://hitran.org/docs/iso-meta/'


def html_formula_to_text(html, for_header=False):
    """Convert HTML formula to text format.

    Data rows:  H<sub>2</sub><sup>16</sup>O  ->  1H2-16O
    Headers:    H<sub>2</sub>O               ->  H2O
    """
    segments = []
    pos = 0
    while pos < len(html):
        if html[pos:pos+5] == '<sup>':
            end = html.index('</sup>', pos)
            segments.append(('sup', html[pos+5:end]))
            pos = end + 6
        elif html[pos:pos+5] == '<sub>':
            end = html.index('</sub>', pos)
            segments.append(('sub', html[pos+5:end]))
            pos = end + 6
        elif html[pos] == '<':
            end = html.index('>', pos)
            pos = end + 1
        else:
            end = html.find('<', pos)
            if end == -1:
                end = len(html)
            segments.append(('text', html[pos:end]))
            pos = end

    result = []
    for i, (seg_type, seg_text) in enumerate(segments):
        if seg_type == 'sub':
            result.append(seg_text)
        elif seg_type == 'sup':
            if seg_text == '+':
                result.append('_p')
            elif i == 0:
                result.append(seg_text)
            else:
                result.append(f'-{seg_text}')
        else:
            result.append(seg_text)

    formula = ''.join(result)

    if not for_header:
        formula = formula.replace('H', '-1H')
        formula = formula.replace('D', '-2H')
        if formula.startswith('-'):
            formula = formula[1:]

    return formula


def convert_scientific(html):
    """Convert HTML scientific notation to E notation.

    9.97317&nbsp;&times;&nbsp;10<sup>-1</sup>  ->  9.97317E-1
    """
    html = html.replace('&nbsp;', ' ').replace('&times;', '\u00d7')
    m = re.match(r'\s*([\d.]+)\s*\u00d7\s*10<sup>([^<]+)</sup>', html)
    if m:
        return f'{m.group(1)}E{m.group(2)}'
    return re.sub(r'<[^>]+>', '', html).strip()


def strip_tags(html):
    """Strip HTML tags and entities, return plain text."""
    text = html.replace('&nbsp;', ' ').replace('&times;', '\u00d7')
    return re.sub(r'<[^>]+>', '', text).strip()


def main():
    print(f'Fetching {URL} ...')
    response = requests.get(URL)
    response.raise_for_status()
    html = response.text

    # Pair each <h4> header with its following <table>
    pattern = r'(<h4>.*?</h4>)|(<table.*?</table>)'
    molecules = []
    current_h4 = None
    for m in re.finditer(pattern, html, re.DOTALL):
        h4, table = m.group(1), m.group(2)
        if h4:
            current_h4 = h4
        elif table and current_h4:
            molecules.append((current_h4, table))
            current_h4 = None

    # Write hitran_meta.txt
    global_ids = set()
    with open('hitran_meta.txt', 'w') as f:
        for h4_html, table_html in molecules:
            h4_inner = re.sub(r'</?h4>', '', h4_html)
            num_str, formula_html = h4_inner.split(':', 1)
            mol_num = int(num_str.strip())
            mol_formula = html_formula_to_text(formula_html.strip(),
                                               for_header=True)

            f.write(f'{mol_num} : {mol_formula}\n')
            f.write('  global  local  Formula  AFGL  Abundance'
                    '  Mass/g.mol^-1  Q(296K)  Q_file  gi\n')
            f.write('  ------  -----  -------  ----  ---------'
                    '  -------------  -------  ------  --\n')

            rows = re.findall(r'<tr>\s*<td.*?</tr>', table_html, re.DOTALL)
            for row_html in rows:
                cells = re.findall(r'<td[^>]*>(.*?)</td>', row_html,
                                   re.DOTALL)
                if len(cells) < 9:
                    continue

                gid = strip_tags(cells[0])
                lid = strip_tags(cells[1])
                formula = html_formula_to_text(cells[2].strip())
                afgl = strip_tags(cells[3])
                abundance = convert_scientific(cells[4])
                mass = strip_tags(cells[5])
                q296 = convert_scientific(cells[6])
                qfile = strip_tags(cells[7])
                gi = strip_tags(cells[8])

                f.write(f'  {gid}  {lid}  {formula}  {afgl}'
                        f'  {abundance}  {mass}  {q296}'
                        f'  {qfile}  {gi}\n')
                global_ids.add(int(gid))

    print(f'Wrote hitran_meta.txt ({len(global_ids)} isotopologues)')

    # Download partition function files
    os.makedirs('Q', exist_ok=True)

    def download_q(n):
        fname = f'q{n}.txt'
        url = f'https://hitran.org/data/Q/{fname}'
        r = requests.head(url, allow_redirects=True)
        if r.status_code != 200:
            return f'URL {url} not found'
        r = requests.get(url)
        r.raise_for_status()
        with open(os.path.join('Q', fname), 'wb') as fout:
            fout.write(r.content)
        return f'Downloaded {fname}'

    with ThreadPoolExecutor(max_workers=10) as executor:
        futures = {executor.submit(download_q, n): n
                   for n in sorted(global_ids)}
        for future in tqdm(as_completed(futures), total=len(futures),
                           desc='Downloading Q files'):
            n = futures[future]
            try:
                future.result()
            except Exception as e:
                print(f'Error downloading q{n}.txt: {e}')


if __name__ == '__main__':
    main()
