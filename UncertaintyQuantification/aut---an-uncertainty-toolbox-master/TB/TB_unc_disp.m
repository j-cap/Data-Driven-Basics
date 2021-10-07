% GenerateDispStr (disp) Testbench
% disp für scalar
clc

% Test Cases:
% TC Sonderfall Null
TC(1).value = unc(0,0.1);
TC(1).expected_str = '0.00(10)';

TC(2).value = unc(0.001351,0.0000001534);
TC(2).expected_str = '0.00135100(15)';

% TC Auffüllen mit Nullern
TC(3).value = unc(1.536e-9,0.0126e-9);
TC(3).expected_str = '1.536e-09(13)';

% TC Negativer Wert
TC(4).value = unc(-144.587,15.1);
TC(4).expected_str = '-145(15)';

% Großer Wert
TC(5).value = unc(22555012,22500);
TC(5).expected_str = '22555000(23000)';

% Wert um 1
TC(6).value = unc(1.1,0.124);
TC(6).expected_str = '1.10(12)';

% Wert um 1
TC(7).value = unc(11.1,1.124);
TC(7).expected_str = '11.1(1.1)';

% Wert kleiner als Unsicherheit und Komma
TC(8).value = unc(0.00001,0.012);
TC(8).expected_str = '0.000(12)';

% TC Sonderfall 1
TC(9).value = unc(1,1);
TC(9).expected_str = '1.0(1.0)';

% Wert kleiner als Unsicherheit und Komma
TC(10).value = unc(55555.0123,612000.012);
TC(10).expected_str = '60000(610000)';

% TC X
TC(11).value = unc(10,1);
TC(11).expected_str = '10.0(1.0)';

% TC X
TC(12).value = unc(1,0.0001);
TC(12).expected_str = '1.00000(10)';

% TC X
TC(13).value = unc(0.1,1.001);
TC(13).expected_str = '0.1(1.0)';

% TC X
TC(14).value = unc(-99.999,12.001);
TC(14).expected_str = '-100(12)';

% TC X
a=unc(0.43,0.001);
TC(15).value = a-a;
TC(15).expected_str = '0(0)';

% TC X
TC(16).value = a/a;
TC(16).expected_str = '1(0)';

% TC X
TC(17).value = unc(0.00001351,0.00000001534);
TC(17).expected_str = '1.3510e-05(15)';

disp ('Start "GenerateDispStr" Test Suite ...');
for Index = 1:1:length(TC)
    expected_str = TC(Index).expected_str;
    DispString = unc.GenerateDispStr(TC(Index).value.value,TC(Index).value.std_unc);

    if ~strcmp(DispString,expected_str)
        disp (DispString);   
        disp (['Testcase: ',num2str(Index),' failed!']);
        disp (['Expected: ',expected_str]);
    end
end
disp ('Test suite finished!');

% disp für Matrix
a = unc([-99.999,100;51,12],[12.001,0.1;8.9,1.2]);
disp ('A 2x2 matrix shall be visible below:');
disp(a);

disp ('A scalar shall be visible below:');
a = unc(3154.457,235.47);
disp(a);

% disp_contibution

a = unc(10,0.1123,'Einen sehr langer Name und ');
b = unc(10,-0.08764e-5,'Kurz');

c=a+b;
disp_contribution(c);

d = unc(10,0.1,'Nur Kurz ');
e = unc(10,0.01e-5,'Kurz');

f=d-e;
g=f+c;
disp_contribution(f);
disp_contribution(g);

msgbox('Please manually verify the proper representation of the results in the command line.','Test description');