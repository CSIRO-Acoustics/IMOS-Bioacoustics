function process_Bluefin_PR(PR_file)

[datetime, pitch, roll]=textread(PR_file,'%s %f %f','delimiter',',','headerlines',1);

hold on
plot(pitch,'b.')
plot(roll,'r.')
legend('pitch','roll')
plot(pitch,'b')
plot(roll,'r')

vector= (sqrt(pitch.^2 + roll.^2));
figure; hold on
plot(vector,'m'); plot(vector,'m.')
