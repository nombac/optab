MODULE string_module

  implicit none
  public

CONTAINS

! https://stackoverflow.com/questions/58938347/how-do-i-replace-a-character-in-the-string-with-another-charater-in-fortran  
    pure recursive function replaceStr(string,search,substitute) result(modifiedString)
        implicit none
        character(len=*), intent(in)  :: string, search, substitute
        character(len=:), allocatable :: modifiedString
        integer                       :: i, stringLen, searchLen
        stringLen = len(string)
        searchLen = len(search)
        if (stringLen==0 .or. searchLen==0) then
            modifiedString = ""
            return
        elseif (stringLen<searchLen) then
            modifiedString = string
            return
        end if
        i = 1
        do
            if (string(i:i+searchLen-1)==search) then
                modifiedString = string(1:i-1) // substitute // replaceStr(string(i+searchLen:stringLen),search,substitute)
                exit
            end if
            if (i+searchLen>stringLen) then
                modifiedString = string
                exit
            end if
            i = i + 1
            cycle
        end do
    end function replaceStr  
  
END MODULE string_module
