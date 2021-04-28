open ss,"<Strain_Stress.dat"; 
@data=<ss>;
close ss;

open dat,">Strain_Stress.dat";
#print"@data";
foreach (@data)
{

   if ($_ =~ m/(lz:\s+\S+\s+Vol:\s+\S+\s+PE:\s+\S+\s+Pxx:\s+\S+\s+Pyy:\s+\S+\s+)(Pzz:\s+)(\S+)(\s+)/x)
   {
       
	    
	    $modify=$3*(-1);
		printf dat "$1$2$modify$4";
	   
   
   }


}