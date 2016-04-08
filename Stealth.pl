############################################################
#
# Stealth.pl
# - Given a transposase sequence, Stealth generates a sequence that is unique 
#   in sequence composition, yet identical in 
#   amino acid composition. Theoretically, the generated 
#   transposase will not be detected by piRNAs. Stealth makes 
#   optimal sequence changes to the provided transposase
#   sequence in terms of tRNA abundance.  Optimal codons  
#   are changed to the next most abundant codon and ties are
#   resolved by changing the input sequence.
#   
#
# by: Patrick Schreiner
# Graduate Student
# Genetics, Genomics, and Bioinformatics
# University of California
# Riverside, CA 92521
# Ph: 847-732-0499
# patrick.schreiner@email.ucr.edu
#
############################################################

#!usr/bin/perl 
use strict;
use File::Basename;
use Getopt::Long;
use Cwd;

my $seqfile;
my $biases;
my $outputname;
my $outputdir;
my $loc;
GetOptions(
	'1:s' => \$seqfile,
	't:s' => \$biases,
	'o:s' => \$outputname,
	'd:s' => \$outputdir,
);

unless ($seqfile and $biases) {
    die ("usage: perl Stealth.pl -1 <Sequence (FASTA format)> -t <Bias Table (Tab-delimited)> -o <OPTIONAL output name> -d <OPTIONAL output directory>\n");
}
unless($outputname)	{
	$loc = index($seqfile,".fasta");
	if($loc == -1)
	{
		$outputname = basename($seqfile,".fa");
	}
	else
	{
		$outputname = basename($seqfile,".fasta");
	}
}
unless($outputdir)	{
	$outputdir=cwd();
	chomp $outputdir;
}
unless(-e $outputname)	{
	my $dir = mkdir($outputname);
	unless($dir)
	{
		die "Failed to create directory $outputdir\n";
	}
}

my $faname = $outputdir . $outputname . '.fa';

my $sequence = parse_fasta($seqfile);
$sequence = uc($sequence);
my @codons=get_codons($sequence);

my($prot,$codes,$values1,$values2,$values3)=&get_biases($biases);
my @aa=@$prot;
my @codon=@$codes;
my @allbias=@$values1;
my @hegbias=@$values2;
my @tebias=@$values3;

chdir $outputname;

my @optimal;
my $final;
my $count=0;
my @stealth;
foreach my $codon(@codons)
{
	open(FILE, ">>Report.txt");
	$count++;
	print FILE "Codon number: $count\n";
	close(FILE);
	@optimal=optimize(codon2aa($codon),\@aa,\@allbias);
	$final=stealth($codon,\@optimal,\@codon,\@allbias);
	push(@stealth,$final);
}

my $output = join("",@stealth);

$faname = $outputname . '.fa';
open(RESULT,">$faname");
print RESULT "> Optimal Stealth Sequence: ";
print RESULT $outputname;
print RESULT "\n";
print RESULT $output;
close(RESULT);
print "\n";

###############################  SUBROUTINES  ###############################

sub parse_fasta	{
	my($file)=@_;
	
	my $i=0;
	my $seq;
	open(FASTA,$file) or die("Could not open file.");
	foreach my $line(<FASTA>)
	{
		if($i == 0)
		{
			$i++;
		}
		else
		{
			$seq.=$line;
		}
	}
	close(FASTA);

	return $seq;
}

sub get_bias	{
	my($cod, $table_cod)=@_;
	
	my @table=@$table_cod;
	my $i=-1;
	my $codon;
	do
	{
		$codon=shift(@table);
		$i++;
	} until ($codon eq $cod);
	return $i;
}

sub get_codons	{
	my($seq)=@_;
	
	for(my $i=0;$i<=(length($seq)-2);$i+=3)
	{
		my $temp=substr($seq,$i,3);
		push(@codons,$temp);
	}
	return @codons;
}

sub get_biases	{
	my($table)=@_;
	
	my @aas;
	my @cods;
	my @biases1;
	my @biases2;
	my @biases3;
	open(BIASTABLE, $table);
	unless( open(BIASTABLE, $table) ) 
	{
		print "Cannot open file: $table\n\n";	
		exit;	
	}
	while(<BIASTABLE>)
	{
		chomp;
		(my $aa, my $cod, my $bias1, my $bias2, my $bias3) = split("\t");
		push(@aas,$aa);
		push(@cods,$cod);
		push(@biases1,$bias1);
		push(@biases2,$bias2);
		push(@biases3,$bias3);
	}
	close(BIASTABLE);
	return (\@aas,\@cods,\@biases1,\@biases2,\@biases3); 
}

sub optimize	{
	my($given_aa,$table_amino,$vals)=@_;
	
	my @table_aa=@$table_amino;
	my @table_values=@$vals;
	my $best=0;
	my $second=0;
	my $test=0;
	for(my $i=0;$i<=$#table_aa;$i+=1)
	{
		if(@table_aa[$i] eq $given_aa)
		{
			$test=@table_values[$i];
			if($test==1)
			{
				$best=$i;
				$second=$i;
			}
			elsif($test>@table_values[$best])
			{
				$second=$best;
				$best=$i;
			}
			elsif($test>@table_values[$second])
			{
				$second=$i;
			}
		} 
	}
	my @positions=($best,$second);
	return @positions; 
}

sub stealth	{
	my($given_cod,$pos,$tab_codon, $vals_interest)=@_;

	my @positions=@$pos;
	my @table_codon=@$tab_codon;
	my @values_of_interest=@$vals_interest;
	my $new_seq;
	my $count=0;
	my $use=0;

	open(FILE, '>>Report.txt');
	if($given_cod eq 'TGA')
	{
		print FILE "No Change: $given_cod\n";
		print FILE "\n";
		$new_seq='TGA';
	}
	elsif($given_cod eq 'TAG')
	{
		print FILE "No Change: $given_cod\n";
		print FILE "\n";
		$new_seq='TAG';
	}
	elsif($given_cod eq 'TAA')
	{
		print FILE "No Change: $given_cod\n";
		print FILE "\n";
		$new_seq='TAA';
	}
	elsif(@positions[0]==@positions[1])
	{
		print FILE "No Change: $given_cod\n";
		print FILE "\n";
		$new_seq=@table_codon[@positions[0]];
	}
	elsif($given_cod eq @table_codon[@positions[0]])
	{
		$new_seq=@table_codon[@positions[1]];
		print FILE "Change from: $given_cod to @table_codon[@positions[1]]\n";
		$use = get_bias($given_cod, \@table_codon);
		$use = @values_of_interest[$use];
		print FILE "Usage value change from: $use to @values_of_interest[@positions[1]]\n";
		print FILE "\n";
	}
	else
	{
		$new_seq=@table_codon[@positions[0]];
		print FILE "Change from: $given_cod to @table_codon[@positions[0]]\n";
		$use = get_bias($given_cod, \@table_codon);
		$use = @values_of_interest[$use];
		print FILE "Usage value change from: $use to @values_of_interest[@positions[0]]\n";
		print FILE "\n";
	}
	close(FILE);
	
	return $new_seq;
}

sub codon2aa {
     my($codon) = @_;
     
	if ( $codon =~ /TCA/i )    { return 'S' }    # Serine
     elsif ( $codon =~ /TCC/i )    { return 'S' }    # Serine
     elsif ( $codon =~ /TCG/i )    { return 'S' }    # Serine
     elsif ( $codon =~ /TCT/i )    { return 'S' }    # Serine
     elsif ( $codon =~ /TTC/i )    { return 'F' }    # Phenylalanine
     elsif ( $codon =~ /TTT/i )    { return 'F' }    # Phenylalanine
     elsif ( $codon =~ /TTA/i )    { return 'L' }    # Leucine
     elsif ( $codon =~ /TTG/i )    { return 'L' }    # Leucine
     elsif ( $codon =~ /TAC/i )    { return 'Y' }    # Tyrosine
     elsif ( $codon =~ /TAT/i )    { return 'Y' }    # Tyrosine
     elsif ( $codon =~ /TAA/i )    { return '_' }    # Stop
     elsif ( $codon =~ /TAG/i )    { return '_' }    # Stop
     elsif ( $codon =~ /TGC/i )    { return 'C' }    # Cysteine
     elsif ( $codon =~ /TGT/i )    { return 'C' }    # Cysteine
     elsif ( $codon =~ /TGA/i )    { return '_' }    # Stop
     elsif ( $codon =~ /TGG/i )    { return 'W' }    # Tryptophan
     elsif ( $codon =~ /CTA/i )    { return 'L' }    # Leucine
     elsif ( $codon =~ /CTC/i )    { return 'L' }    # Leucine
     elsif ( $codon =~ /CTG/i )    { return 'L' }    # Leucine
     elsif ( $codon =~ /CTT/i )    { return 'L' }    # Leucine
     elsif ( $codon =~ /CCA/i )    { return 'P' }    # Proline
     elsif ( $codon =~ /CCC/i )    { return 'P' }    # Proline
     elsif ( $codon =~ /CCG/i )    { return 'P' }    # Proline
     elsif ( $codon =~ /CCT/i )    { return 'P' }    # Proline
     elsif ( $codon =~ /CAC/i )    { return 'H' }    # Histidine
     elsif ( $codon =~ /CAT/i )    { return 'H' }    # Histidine
     elsif ( $codon =~ /CAA/i )    { return 'Q' }    # Glutamine
     elsif ( $codon =~ /CAG/i )    { return 'Q' }    # Glutamine
     elsif ( $codon =~ /CGA/i )    { return 'R' }    # Arginine
     elsif ( $codon =~ /CGC/i )    { return 'R' }    # Arginine
     elsif ( $codon =~ /CGG/i )    { return 'R' }    # Arginine
     elsif ( $codon =~ /CGT/i )    { return 'R' }    # Arginine
     elsif ( $codon =~ /ATA/i )    { return 'I' }    # Isoleucine
     elsif ( $codon =~ /ATC/i )    { return 'I' }    # Isoleucine
     elsif ( $codon =~ /ATT/i )    { return 'I' }    # Isoleucine
     elsif ( $codon =~ /ATG/i )    { return 'M' }    # Methionine
     elsif ( $codon =~ /ACA/i )    { return 'T' }    # Threonine
     elsif ( $codon =~ /ACC/i )    { return 'T' }    # Threonine
     elsif ( $codon =~ /ACG/i )    { return 'T' }    # Threonine
     elsif ( $codon =~ /ACT/i )    { return 'T' }    # Threonine
     elsif ( $codon =~ /AAC/i )    { return 'N' }    # Asparagine
     elsif ( $codon =~ /AAT/i )    { return 'N' }    # Asparagine
     elsif ( $codon =~ /AAA/i )    { return 'K' }    # Lysine
     elsif ( $codon =~ /AAG/i )    { return 'K' }    # Lysine
     elsif ( $codon =~ /AGC/i )    { return 'S' }    # Serine
     elsif ( $codon =~ /AGT/i )    { return 'S' }    # Serine
     elsif ( $codon =~ /AGA/i )    { return 'R' }    # Arginine
     elsif ( $codon =~ /AGG/i )    { return 'R' }    # Arginine
     elsif ( $codon =~ /GTA/i )    { return 'V' }    # Valine
     elsif ( $codon =~ /GTC/i )    { return 'V' }    # Valine
     elsif ( $codon =~ /GTG/i )    { return 'V' }    # Valine
     elsif ( $codon =~ /GTT/i )    { return 'V' }    # Valine
     elsif ( $codon =~ /GCA/i )    { return 'A' }    # Alanine
     elsif ( $codon =~ /GCC/i )    { return 'A' }    # Alanine
     elsif ( $codon =~ /GCG/i )    { return 'A' }    # Alanine
     elsif ( $codon =~ /GCT/i )    { return 'A' }    # Alanine
     elsif ( $codon =~ /GAC/i )    { return 'D' }    # Aspartic Acid
     elsif ( $codon =~ /GAT/i )    { return 'D' }    # Aspartic Acid
     elsif ( $codon =~ /GAA/i )    { return 'E' }    # Glutamic Acid
     elsif ( $codon =~ /GAG/i )    { return 'E' }    # Glutamic Acid
     elsif ( $codon =~ /GGA/i )    { return 'G' }    # Glycine
     elsif ( $codon =~ /GGC/i )    { return 'G' }    # Glycine
     elsif ( $codon =~ /GGG/i )    { return 'G' }    # Glycine
     elsif ( $codon =~ /GGT/i )    { return 'G' }    # Glycine
     else {
	print STDERR "Invalid codon \"$codon\"!!\n";
	exit;
     }
}