# !/bin/perl

my $NEXEC = 2;
my $MNAME = "./result/accessi_m3.dat";
my $FNAME = "./result/accessi_f3.dat";
my $nMNAME;
my $nFNAME;
my $count = 1;

while ($count <= 7) {

    $nMNAME = "./result/accessi_m3_".$count.".dat";
    $nFNAME = "./result/accessi_f3_".$count.".dat";
    
    print $count."\n";
    
    system("./accessi3.exe >& esito & wait");
    system("mv $MNAME $nMNAME");
    system("mv $FNAME $nFNAME");

    $count++;
}
