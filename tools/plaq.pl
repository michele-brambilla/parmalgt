# !/bin/perl

# Primi valori della placchetta. Da invocare con il numero di colori
# e il numero di flavour. Nel caso Nf != 0 restituiece solo i primi 
# due coefficienti della placchetta (qualcosa non torna sul terzo..)

$Nc=-1;
$Nf=-1;

($Nc, $Nf) = @ARGV;

$len=@ARGV;

if($len!=2){
print "\nERRORE\n";
print "La sintassi del comando e' plaq.pl <Nc> <Nf>\n\n";
exit;
}

print "\n";
print " Calcola i primi valori della serie perturbativa della placchetta\n";
print "\n";

print "#######################\n";
print "\t"."Nc = ".$Nc."\n";
print "#######################\n";

my $c1g = ($Nc*$Nc-1)/(8*$Nc);
my $c2g = ($Nc*$Nc-1)*(0.0051069297-1.0/(128.0*$Nc*$Nc));
my $c3g = ($Nc*$Nc-1)*(0.0023152583/($Nc*$Nc*$Nc)-0.002265487/$Nc+0.000794223*$Nc);

my $h2  = -0.6929202;
my $h30 =  0.01867;
my $h31 =  0.1479;
my $h32 =  0.04159747;

$h2  = $h2/1000.0;
$h30 = $h30/10000.0;
$h31 = $h31/10000.0;
$h32 = $h32/10000.0;

#print $h30."\n";
#print $h31."\n";
#print $h32."\n";

my $c1f = 0;
my $c2f = ($Nc*$Nc-1)*$h2*$Nf/$Nc;
my $c3f = ($Nc*$Nc-1)*($h30*$Nf+$h31*$Nf/($Nc*$Nc)+$h32*$Nf*$Nf/$Nc);

#print $c1f."\n";
#print $c2f."\n";
#print $c3f."\n";

#print $c1g+$c1f."\n";
#print $c2g+$c2f."\n";
#print $c3g+$c3f."\n";

my $c1 =                     (2.0*$Nc)*($c1g + $c1f);
my $c2 =           (2.0*$Nc)*(2.0*$Nc)*($c2g + $c2f);
my $c3 = (2.0*$Nc)*(2.0*$Nc)*(2.0*$Nc)*($c3g + $c3f);

print "c1 = ".$c1."\n";
print "c2 = ".$c2."\n";
#if($Nf==0)
#   {
     print "c3 = ".$c3."\n";
#   }
