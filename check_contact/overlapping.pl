use warnings;
use strict;
use Parallel::ForkManager;

my $forkNo = `lscpu|grep "^CPU(s)"|awk '{print \$2}'`;
print "forkNo: $forkNo\n";
my $pm = Parallel::ForkManager->new("$forkNo");

my $datafile = "box_relax_done.data";
open my $file_content,"< $datafile" or die("Can't open $datafile");
my @content = <$file_content>;
close ($file_content);
my @dataobj = grep {if(m/^(\d+)\s+(\d+)\s+(\d+)\s+(\-*\d\.?\d*)\s+(\-?\d+\.\d+)\s+(\-?\d+\.\d+)\s+(\-?\d+\.\d+)\s+\-?\d\s+\-?\d\s+\-?\d$/){$_ = [$1,$2,$3,$4,$5,$6,$7];}} @content;
if(! @dataobj) {die "No coorobj information\n";}

open my $result,"> before_MD.txt";
for my $atom1(0..$#dataobj)
	{
		print "doing atom $atom1\n";
		my $loop2=$atom1+1;
$pm->start and next;
		for my $atom2($loop2..$#dataobj)
		{
		my 	$firstx = $dataobj[$atom1][4];
		my 	$firsty = $dataobj[$atom1][5];
		my 	$firstz = $dataobj[$atom1][6];
		my  $secondx = $dataobj[$atom2][4];
		my  $secondy = $dataobj[$atom2][5];
		my  $secondz = $dataobj[$atom2][6];
		my $distance = ((($secondx-$firstx)*($secondx-$firstx))+(($secondy-$firsty)*($secondy-$firsty))+(($secondz-$firstz))*($secondz-$firstz));
		#print "$dataobj[$atom1][0] & $dataobj[$atom2][0] distance == $distance\n";
		if ($distance < 0.04)
		{
			print "$dataobj[$atom1][0] $dataobj[$atom1][1] $dataobj[$atom1][2] $dataobj[$atom1][3] & $dataobj[$atom2][0] $dataobj[$atom2][1] $dataobj[$atom2][2] $dataobj[$atom2][3] are overlapping (Distance == $distance)\n";
		}
		}
$pm-> finish;
	}
$pm->wait_all_children;
close ($result);
