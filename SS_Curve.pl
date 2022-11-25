use warnings;
use strict;

#here-doc for perl
my $pyfile = <<'END_MESSAGE';
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
savename = "Stress-Strain.jpg"
data = np.genfromtxt("./SS.txt", names=True)
plt.gca().xaxis.set_major_locator(MaxNLocator(integer=True))
plt.plot(data['strain'], data['stress'], label='stress-strain curve',marker='o',markevery= 5, linewidth=0.25, markersize=1)
plt.legend()
plt.xlabel('Strain')
plt.ylabel('Stress (GPa)')
plt.savefig(savename)
plt.show()
END_MESSAGE

open my $ss,"> ./plot.py";
printf $ss "$pyfile\n";
close $ss; 

if( -e "./HEA_Tension/Strain_Stress.dat"){
    my @strain = `cat ./HEA_Tension/Strain_Stress.dat| awk '{print \$2}'`;
    my @stress = `cat ./HEA_Tension/Strain_Stress.dat| awk '{print \$12}'`;
    chomp (@strain);
    chomp (@stress);
    my $refL = $strain[0];
    `touch SS.txt`;
    `echo 'strain stress' > SS.txt`;
     for (0..$#strain){
        $strain[$_] = ($strain[$_] - $refL)/$refL;
       # print "check: $_ $strain[$_] $stress[$_]\n";
        `echo '$strain[$_] $stress[$_]' >> SS.txt`;
    }    
    `python ./plot.py`;  
}
else{
    die "No Strain_Stress.dat for plotting stress-strain curve.\n";
}
