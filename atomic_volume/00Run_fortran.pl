=b
convert f77 to f90 or f95
=cut
use warnings;
use strict;
use Data::Dumper;
use Cwd;
use POSIX;

system("gfortran -o Getvol.x Getvol.f");
die "compiling failed" if($?);

`rm -f input.cfg`;#remove old input cfg

my @cfg_files = <*.cfg>;
`rm -f output.dat`;
`touch output.dat`;
for (@cfg_files){
    chomp;
    print "Getting the system volume of $_\n";
    `rm -f input.cfg`;
    `cp $_ input.cfg`;  
    sleep(1);
    system("./Getvol.x");
    system("cat ./cfg_vol.dat >> output.dat");
} 
