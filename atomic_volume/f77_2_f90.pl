=b
convert f77 to f90 or f95
dos2unix:
perl -p -e 's/\r$//' print.sh > print_unix.sh
=cut
use warnings;
use strict;
use Data::Dumper;
use Cwd;
use POSIX;
my @fortran_files = <*.f>;
chomp @fortran_files;
my $out_folder = "F90";
`rm -rf ./$out_folder`;
`mkdir -p ./$out_folder`;
for (@fortran_files){
    $_ =~ /(.+)\.f/;
    $1  =~ s/^\s+|\s+$//;
    my $prefix = $1;
    print "\$prefix: $prefix\n";
    `touch ./$out_folder/$prefix.f90`;
    my @temp = `cat $_`;
    chomp @temp;    
    for my $ln (@temp){
        $ln =~ s/^\s+|\s+$//;
        `echo "$ln" >> ./$out_folder/$prefix.f90`;
    }
}