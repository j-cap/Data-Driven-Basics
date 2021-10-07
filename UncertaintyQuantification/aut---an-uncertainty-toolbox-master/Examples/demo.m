% basic examples for uncertainty calculations

a=unc(100,1);
disp('Create variables a and b with estimate 100 and standard uncertainty of 1:')
disp('a=unc(100,1)')
a=unc(100,1)    
disp('b=unc(100,1)');
b=unc(100,1)   

disp('The value in parenthesis represents the standard uncertainty');
disp(['in terms of the least significant digits of the variable!',13])
pause;
disp('Now let us calculate the difference c=a-b:')

c=a-b           % difference is zero, but has uncertainty
disp(['The difference is zero, but has uncertainty of 1.41!',13])

pause;
disp('Now let us calculate the ratio r=a/b:')

r=a/b           % ratio is one, but has uncertainty
disp(['The ratio is one, but has uncertainty of 0.0141!'])
disp('NOTE: The value in parenthesis represents the standard uncertainty');
disp(['in terms of the least significant digits of the variable!',13])

pause;
disp('Now assign that the two variables a and b are equivalent')
disp('b=a')
b=a

pause;
disp('and let us calculate the difference c=a-b:')

c=a-b           % difference is zero and has zero uncertainty
disp('Now, the difference is zero and has zero uncertainty!')


pause;
disp('Now let us calculate the ratio r=a/b:')

r=a/b           % ratio is one, and has zero uncertainty
disp(['The ratio is one and has uncertainty of 0!',13])

pause;
disp('For further examples have a look into Bridge.m and Hall.m!')
