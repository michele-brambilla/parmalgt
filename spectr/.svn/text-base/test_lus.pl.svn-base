#!/usr/bin/perl

#definizioni:
$NUM=10;                                   #numero di iterazioni per ogni configurazione
$NCONF=20;                                 #numero di configurazioni
$CONFNM="w_conf_6x6x6x6_a0.0_t0.01_o4_";   #nome configurazioni
$CONFPATH="./";                            #percorso relativo

open(INFILE,"testconf_1.txt") || die "File di configurazione non esistente";
@linee = <INFILE>;
close(INFILE);

for ($i = 1; $i <= $NCONF; $i++)
{
    $nome = "conf\t".$CONFPATH.$CONFNM.$i.".dat";
    $linee[4] = $nome."\n";
    for($j = 0; $j < $NUM; $j++){
	$rnd = "seed\t".int(rand 99999999);
	$linee[3] = $rnd."\n";
	open(OUFILE,">testconf.txt");
	print OUFILE @linee;
	close(OUFILE);
	system("./test_lus.exe >> ./esito_lus.txt 2>> ./esito_lus.log");
	print "Eseguito ", 100*(($i-1)*$NUM+$j)/($NUM*$NCONF), "%\n";
    }
}
