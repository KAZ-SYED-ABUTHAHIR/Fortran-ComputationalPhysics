PROGRAM String
    implicit none
    character,dimension(5) :: vowels = (/ 'a','e','i','o','u' /)
    character (len=12) :: string2 = "Hello Again !"
    character(19) :: string3 = " What is your name?"
    character(30) :: string4
    string4 = string2(1:5)//string3
    print*, string4
END PROGRAM String