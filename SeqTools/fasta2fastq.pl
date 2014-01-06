#!/usr/bin/perl

=head1 NAME

 fasta2fastq.pl
 A tools to change fasta to fastq format.

=cut

=head1 SYPNOSIS

 fasta2fastq.pl [-h] -f <fastafile> -q <qualityfile> 
                    [-t <ngs_technology>] [-d <default_qscore>]
                    [-S <seq_number>]


=head2 I<Flags:>

=over

=item -f

B<fastafile>              sequence file in fasta format (mandatory)

=item -q

B<qualityfile>            qscore file in qual format (optional)

=item -t

B<ngs_technology>         ngs_technology used (sanger by default)

=item -d

B<default_qscore>         value to use as qscore when none is used (15 default)

=item -s

B<seq_number>             slice the fasta and quality files before run the
                          script (charge less sequences in the computer memory)
=item -Q

B<quiet>                  run and dont print any message.

=item -V

B<version>                print script version

=item -h

B<help>                   print the help

=back

=cut

=head1 DESCRIPTION

 This script change fasta format to fastq, using the qscores supplied by the 
 -q option or by -d option (as default value).

 Different NGS technologies can be used to modify the coding of the qscore into
 the fastq files:

 + sanger, qscore conversion to ascii (qscore + 33)
 + 454, qscore conversion to ascii (qscore + 33)
 + solexa1.0/illumina1.0, qscore conversion to ascii (qscore + 54)
 + solexa1.3/illumina1.3, qscore conversion to ascii (qscore + 64)
 + solexa1.5/illumina1.5, qscore conversion to ascii (qscore + 64*)
   (* values 0,1 and 2 are not qscores, they coding the following messages:
   0 and 1 = unused, 2 = Read Segment Quality Control)
   Use solexa1.3/illumina1.3 for this case.

 Output will have the name of the input file with the new extension .fastq

=cut

=head1 AUTHORS

  Aureliano Bombarely Gomez.
  (ab782@cornell.edu).

=cut

=head1 METHODS

 fasta2fastq.pl

=cut

use strict;
use warnings;
use autodie;

use Getopt::Std;
use Bio::SeqIO;
use Bio::Seq::PrimaryQual;
use File::Temp qw/ tempfile /;

our $VERSION = '0.03';
$VERSION = eval $VERSION;

our ($opt_f, $opt_q, $opt_t, $opt_d, $opt_s, $opt_Q, $opt_V, $opt_h);
getopts("f:q:t:d:s:QVh");
if (!$opt_f && !$opt_q && !$opt_t && !$opt_d && !$opt_s && !$opt_Q && !$opt_V) {
    print "There are n\'t any tags. Print help\n\n";
    help();
}
if ($opt_h) {
    help();
}
if ($opt_V) {
    print STDERR "\nfasta2fastq.pl VERSION: $VERSION\n\n";
    exit();
}

## Get the arguments and check them

my $fasta = $opt_f || 
    die("INPUT ARGS ERROR: Input option was not supplied (-f <fasta_file>).\n");

my $qual = $opt_q;
my $tech = $opt_t || 'sanger';
my $defqual = $opt_d || '15';
my $slices = $opt_s;

my %tech = ( 
    'sanger'      => 33,
    '454'         => 33,
    'solexa1.0'   => 54,
    'illumina1.0' => 54,
    'solexa1.3'   => 64,
    'illumina1.3' => 64,
    );

my $date = `date`;
chomp($date);

print STDERR "\n************************************************************\n";
print STDERR "* fasta2fastq script begins ($date) *\n";
print STDERR "************************************************************\n";

## Now it will extract the fasta sequences and put them into a hash


my $outfile = $fasta;
if ($outfile =~ m/\.\w+$/) {
    $outfile =~ s/\.\w+$/\.fastq/;
}
else {
    $outfile .= '.fastq';
}


if (defined $slices) {

    ## Split the files

    my @tempfasta = ();
    my @tempqual = ();

    my $fastaio = Bio::SeqIO->new( -file => $fasta, -format => 'fasta' );

    my $newfastaio;
    my $sq = 0;

    unless ($opt_Q) {
	print STDERR "FILE SLICING OPTION IS ENABLED ($slices)\n\n";
    }

    while( my $seq = $fastaio->next_seq()) {
    
	unless (defined $newfastaio) {
	    my ($tempfh, $tempfile) = tempfile("$fasta"."_XXXXXX", 
					       SUFFIX => '.fasta',
					       UNLINK => 1,
					       TMPDIR => 1,
		);
	    close($tempfh);

	    $newfastaio = Bio::SeqIO->new(-file => ">$tempfile", 
					  -format => 'fasta' );

	    push @tempfasta, $tempfile;
	}
	$newfastaio->write_seq($seq);
	$sq++;

	if ($sq == $slices) {
	    undef($newfastaio); 
	    $sq = 0;
	}
    }

    if (defined $qual) {
	my $qualio = Bio::SeqIO->new( -file => $qual, -format => 'qual' );

	my $newqualio;
	my $qq = 0;

	while( my $qua = $qualio->next_seq()) {
    
	    unless (defined $newqualio) {
		my ($tempfh, $tempfile) = tempfile("$qual"."_XXXXXX", 
						   SUFFIX => '.qual',
						   UNLINK => 1,
						   TMPDIR => 1,
		    );
		close($tempfh);

		$newqualio = Bio::SeqIO->new(-file => ">$tempfile", 
					     -format => 'qual' );

		push @tempqual, $tempfile;
	    }
	    $newqualio->write_seq($qua);
	    $qq++;

	    if ($qq == $slices) {
		undef($newqualio); 
		$qq = 0;
	    }
	} 
    }

    unless ($opt_Q) {
	my $count = scalar(@tempfasta);
	print STDERR "$count FILES HAVE BEEN CREATED\n\n";
    }

    my $index = 0;
    my @outfiles = ();
    
    foreach my $sl_fasta (@tempfasta) {
    
	my $sl_qual = $tempqual[$index];

	my %seq = get_fasta($sl_fasta);

	my %qual = get_qual(\%seq, $sl_qual, $defqual);

	my $sl_outfile = $sl_fasta;
	if ($sl_outfile =~ m/\.\w+$/) {
	    $sl_outfile =~ s/\.\w+$/\.fastq/;
	}
	else {
	    $sl_outfile .= '.fastq';
	}
	
	print_fastq(\%seq, \%qual, $tech, $sl_outfile);
	$index++;
	push @outfiles, $sl_outfile;
    }
    my $catfiles = join(' ', @outfiles);
    
    system("cat $catfiles > $outfile");

    unless ($opt_Q) {
	print STDERR "JOINING OUTPUT FILES INTO $outfile\n\n";
    }
}
else {

    my %seq = get_fasta($fasta);

    my %qual = get_qual(\%seq, $qual, $defqual);
    
    print_fastq(\%seq, \%qual, $tech, $outfile);
}

my $date2 = `date`;
chomp($date2);

print STDERR "\n**********************************************************\n";
print STDERR "* fasta2fastq script ends ($date2) *\n";
print STDERR "**********************************************************\n";


=head2 help

  Usage: help()
 
  Desc: print help of this script
 
  Ret: none
 
  Args: none
 
  Side_Effects: exit of the script
 
  Example: if (!@ARGV) {
               help();
           }

=cut

sub help {
  print STDERR <<EOF;
  $0:

    Description:

    This script change fasta format to fastq, using the qscores supplied by the 
    -q option or by -d option (as default value).

    Different NGS technologies can be used to modify the coding of the qscore 
    into the fastq files:

     + sanger, qscore conversion to ascii (qscore + 33)
     + 454, qscore conversion to ascii (qscore + 33)
     + solexa1.0/illumina1.0, qscore conversion to ascii (qscore + 54)
     + solexa1.3/illumina1.3, qscore conversion to ascii (qscore + 64)
     + solexa1.5/illumina1.5, qscore conversion to ascii (qscore + 64*)
       (* values 0,1 and 2 are not qscores, they coding the following messages:
       0 and 1 = unused, 2 = Read Segment Quality Control)
       Use solexa1.3/illumina1.3 for this case.

     Output will have the name of the input file with the new extension .fastq

    Usage:
     
      fasta2fastq.pl [-h] -f <fastafile> -q <qualityfile> 
                    [-t <ngs_technology>] [-d <default_qscore>]

    Examples:

      fasta2fastq.pl -f myfile.fasta -q myfile.qual 

    Flags:

     -f <fastafile>              sequence file in fasta format (mandatory)
     -q <qualityfile>            qscore file in qual format (optional)
     -t <ngs_technology>         ngs_technology used (sanger by default)
     -d <default_qscore>         value to use as qscore when none is used 
                                 (15 default)
     -S <seqnumber>              slice the fasta file and qual file before
                                 parse (reduce the memory requeriments)
     -V <version>                print version
     -h <help>                   print the help
    

EOF
exit (1);
}

=head2 get_fasta

  Usage: my %fasta = get_fasta($fastafile);
 
  Desc: Parse the fasta file and return a hash with key=id, value=Bio::Seq
 
  Ret: A hash with key=id and value=Bio::Seq
 
  Args: A fasta filename
 
  Side_Effects: print status if -Q is used.
 
  Example: my %fasta = get_fasta($fastafile);

=cut

sub get_fasta {
    my $fasta = shift;
    
    my %seq = ();

    my $s = 0;
    
    unless ($opt_Q) {
	print STDERR "\n\tPARSING FASTA FILE $fasta\n\n";
    }

    my $fastaio = Bio::SeqIO->new( -file => $fasta, -format => 'fasta');
    
    while (my $seqobj = $fastaio->next_seq()) {
	my $id = $seqobj->display_id();
	$s++;
	
	unless ($opt_Q) {
	    print STDERR "\t\tLoading sequence $s with ID: $id     \r";
	}

	$seq{$id} = $seqobj;
    }

    unless ($opt_Q) {
	print STDERR "\n";
    }

    return %seq;
}

=head2 get_qual

  Usage: my %qual = get_fasta($seqhref, $qualfile, $defqual);
 
  Desc: Parse the fasta file and return a hash with key=id, 
        value=Bio::Seq::PrimaryQual
 
  Ret: A hash with key=id and value=Bio::Seq::PrimaryQual
 
  Args: A hashref with key=ID and value=Bio::Seq
        A qual filename or a default quality value
 
  Side_Effects: print status if -Q is used.
 
  Example: my %qual = get_qual(\%seqs, $fastafile);
           my %qual = get_qual(\%seqs, undef, $defqual);

=cut

sub get_qual {
    my $seqhref = shift;
    my $qual = shift;
    my $defqual = shift;

    my %seq = %{$seqhref};

    my %qual = ();
    my $q = 0;

    if (defined $qual) {
	
	unless ($opt_Q) {
	    print STDERR "\n\tPARSING QUAL FILE $qual\n\n";
	}
	
	my $qualio = Bio::SeqIO->new( -file => $qual, -format => 'qual');
    
	while (my $qualobj = $qualio->next_seq()) {
	    $q++;
	    my $id = $qualobj->display_id();
	 
	    unless ($opt_Q) {
		print STDERR "\t\tLoading quality $q with ID: $id     \r";
	    }
	    $qual{$id} = $qualobj;
	}
    }
    else {
	
	unless ($opt_Q) {
	    print STDERR "\n\tNO QUAL FILE WAS USED. DEFAULT QUAL VALUES\n\n";
	}
	
	foreach my $id (keys %seq) {
	    $q++;
	    my @qual = ();
	    my @seqnt = split(//, $seq{$id}->seq());
	    foreach my $nt (@seqnt) {
		push @qual, $defqual;
	    }
	    unless ($opt_Q) {
		print STDERR "\t\tCreating quality $q for ID: $id     \r";
	    }
	    
	    my $qualobj = Bio::Seq::PrimaryQual->new( -qual => join(' ', @qual),
						      -id   => $id,
		);
	    $qual{$id} = $qualobj;
	}
    }
    unless ($opt_Q) {
	print STDERR "\n";
    }

    return %qual;
}


=head2 print_fastq

  Usage: print_fastq(\%fasta, \%qual, $tech, $outfile);
 
  Desc: Print a fastq file
 
  Ret: None
 
  Args: \%fasta, a hash ref. with key=ID and value=Bio::Seq
        \%qual, a hash re. with key=ID and value=Bio::Seq::PrimaryQual
        $tech, a technology type used.
        $outfile, an outfilename
 
  Side_Effects: print status if -Q is used.
 
  Example: print_fastq(\%fasta, \%qual, $outfile);

=cut

sub print_fastq {
    my $seq_href = shift;
    my $qual_href = shift;
    my $tech = shift;
    my $outfile = shift;

    my %seq = %{$seq_href};
    my %qual = %{$qual_href};

    unless ($opt_Q) {
	print STDERR "\n\tMERGING SEQ-QUAL.PRINTING FASTQ FILE ($outfile)\n\n";	
    }

    my $fastqio = Bio::SeqIO->new( -file => ">$outfile", -format => 'fastq');
    if ($tech =~ m/(illumina|solexa)/) {

	if ($tech =~ m/1\.0$/) {
	    $fastqio->variant('solexa');
	}
	else {
	    $fastqio->variant('illumina');
	}
    }
    else {
	$fastqio->variant('sanger');
    }
    
    my $f = 0;
    foreach my $id (sort keys %seq) {
	
	$f++;
	my $qual = Bio::Seq::Quality->new( -id   => $id,
					   -seq  => $seq{$id}->seq(),
					   -qual => $qual{$id}->qual(),
	    );

	unless ($opt_Q) {
	    print STDERR "\t\tPrinting fastq $f for ID: $id     \r";
	}
	
	$fastqio->write_fastq($qual);
    }
    unless ($opt_Q) {
	print STDERR "\n\n";
    }
}
