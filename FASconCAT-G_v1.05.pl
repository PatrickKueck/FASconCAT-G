#!/usr/bin/perl
use strict;
use File::Copy;
use Tie::File;
use Term::ANSIColor qw(:constants);
use Getopt::Std;

######################################### ENTER THE PROTTEST SOFTWARE NAME OF YOUR RUNNING SYSTEM IN SINGLE QUOTATION MARKS ''
# DEFINEMENT OF ACTUAL PROTTEST SOFTWARE NAME
my $prottest = 'prottest-3.3.jar' ;
##############################


# written by patrick kück, zentrales forschungsmuseum alexander koenig, bonn, germany
# email: patrick_kueck@web.de

# updated on 12th,january	, 2009 by patrick kueck
# updated on 25th may		, 2009 by patrick kueck
# updated on 28th may		, 2009 by patrick kueck
# updated on  7th june		, 2009 by patrick kueck
# updated on  8th june		, 2009 by patrick kueck
# updated on 28th september	, 2009 by patrick kueck
# updated on 16th october	, 2009 by patrick kueck
# updated on 17th november	, 2009 by patrick kueck
# updated on  7th january	, 2010 by patrick kueck
# updated on  7th april		, 2010 by patrick kueck
# updated on  8th april		, 2010 by patrick kueck -> FcC v1.0
# updated on      april		, 2014 by patrick kueck -> FcC-G v1.0
# updated on 30th september	, 2014 by patrick kueck -> FcC-G v1.01 Implementation of additional Screen Output (%Missing Data and Single Infile Sequence Ranges in the concatenated supermatrix) as well as some output result file format changes
# updated on 16th november	, 2014 by patrick kueck ->input file endings line 1529-1531->FcC v1.03
# updated on 22nd november	, 2016 by patrick kueck -> concatenation and print out of supermatrix ranges due to sorted infile names and not randomly. Debuged read in error of .FASTA or .fasta named files
# updated on 22nd march     , 2021 by patrick kueck -> FcC-G v1.05 Including an additional character coding option for absent sequence data (-g). With -g, absent sequence data is filled by '-' rather than 'N/X' during the concatenation process
####################################################### START #######################################################################################



my %outfile_name = (
						'supermatrixFAS'	=> 'FcC_smatrix.fas',		#	Supermatrix outfile in FASTA format
						'supermatrixPHY'	=> 'FcC_smatrix.phy',		#	Supermatrix outfile in PHYLIP format
						'supermatrixNEX'	=> 'FcC_smatrix.nex',		#	Supermatrix outfile in NEXUS format
						'structure'			=> 'FcC_structure.txt',		#	Structure sequence info file
						'info'				=> 'FcC_info.xls',			#	Sequence states & concatenation info file
						'prottest'			=> $prottest,				#	Name of defined prottestversion
) ;


my @parameter_all	= ( 'YES',			'NO'								) ;	# Parameter option of data concatenation					-> 'YES' -> concatenate all files;	'NO' -> concatenate defined files
my @parameter_info	= ( 'YES',			'NO'								) ;	# Parameter option of additional sequence state analysis	-> 'YES' -> perform and print ;		'NO' -> print only basal sequence info
my @parameter_tra	= ( 'NO',			'NUC to AA','AA to NUC'				) ;	# Parameter option of sequence transaltion					-> 'NUC to AA' -> translate nucleotide data to amino acid data; 'AA to NUC' -> translate amino acid data to nucleotide data
my @parameter_phy	= ( 'NO',			'STRICT',	'RELAXED'				) ;	# Print supermatrix in PHYLIP format
my @parameter_nex	= ( 'NO',			'BLOCK',	'MrBAYES'				) ;	# Print supermatrix in NEXUS format
my @parameter_fas	= ( 'YES',			'NO'								) ;	# Print supermatrix in FASTA format
my @parameter_con	= ( 'NO',			'Freq',		'Maj',		'Strict'	) ;	# Perform consensus sequence of defined sequence blocks
my @parameter_file	= ( 'Supermatrix',	'Convert',	'Supermatrix/Convert'	) ; # Concatenate, Convert, or both define by associated option	-> 'Supermatrix': concatenation; Convert: Converting single files without adding missing sequences; 'Supermatrix/Convert': Concatenated single files and convert single files with included missing sequence taxa
my @parameter_3rd	= ( 'Remain',		'Reject'							) ; # Parameter option of 3rd sequence position exclusion		-> 'Remain' -> remove third position;	'Reject' -> keep third position
my @parameter_ryc	= ( 'NO',			'All',		'3rd'					) ; # Parameter option of RY coding								-> 'NO' -> no RY coding;	'All' -> RY coding of complete sequence;	'3rd' -> RY coding of 3rd sequence positions
my @parameter_par	= ( 'NO',			'YES'								) ;	# Print parsimonious sites as extra msa file				-> 'NO' -> no print;	'YES' -> print
my @parameter_ren	= ( 'NO',			'YES'								) ;	# Rename taxon names of given infiles						-> 'NO' ;	'YES'
my @parameter_prt	= ( 'NO',			'Supermatrix'						) ;	# Print OUT additional partition files if concatenated		-> 'NO' ;	'Supermatrix' means yes for concatenatde files
my @parameter_pro	= ( 'NO',			'Supermatrix'						) ;	# Start prottest analyses für aa data						-> 'NO' ;	'YES'
my @parameter_mis	= ( 'Missing',		'Indel'								) ;	# Replace missing sequences with							-> 'Missing': 'X' or 'N' ; 'Indel': '-'

##############################
## READ IN ARGV and assign user defined parameter options
&argv_handling (
					\@parameter_all,	# List of parameter options of data concatenation					-> IN (defined) / OUT (changed)
					\@parameter_info,	# List of parameter options of info print out						-> IN (defined) / OUT (changed)
					\@parameter_phy,	# List of supermatrix parameter options of PHYLIP print OUT			-> IN (defined) / OUT (changed)
					\@parameter_nex,	# List of supermatrix parameter options of NEXUS print OUT			-> IN (defined) / OUT (changed)
					\@parameter_con,	# List of consensus options											-> IN (defined) / OUT (changed)
					\@parameter_file,	# List of filehandling options										-> IN (defined) / OUT (changed)
					\@parameter_fas,	# List of supermatrix parameter options of FASTA print OUT			-> IN (defined) / OUT (changed)
					\@parameter_3rd,	# List of third position handling									-> IN (defined) / OUT (changed)
					\@parameter_ryc,	# List of RY coding													-> IN (defined) / OUT (changed)
					\@parameter_tra,	# List of Sequence transaltion options								-> IN (defined) / OUT (changed)
					\@parameter_par,	# Print parsimonious sites as extra msa file						-> IN (defined) / OUT (changed)
					\@parameter_ren,	# Rename taxon names of given ifiles								-> IN (defined) / OUT (changed)
					\@parameter_prt,	# Print partition files for concatenated data						-> IN (defined) / OUT (changed)
					\@parameter_pro,	# Start prottest analyses für aa data								-> IN (defined) / OUT (changed)
					\@parameter_mis,	# Replace missing gene sequences by '-' instead using 'X' or 'N'	-> IN (defined) / OUT (changed)
					\%outfile_name,		# Outfilename of output option										-> IN (defined) / OUT (unchanged)
) ;
##############################



##############################
## Open MENU
&menu ;
##############################



##############################
## Open Parameter setup
&parameter(
			\@parameter_all,	# List of parameter options of data concatenation					-> IN (defined) / OUT (changed)
			\@parameter_info,	# List of parameter options of info print out						-> IN (defined) / OUT (changed)
			\@parameter_phy,	# List of supermatrix parameter options of PHYLIP print OUT			-> IN (defined) / OUT (changed)
			\@parameter_nex,	# List of supermatrix parameter options of NEXUS print OUT			-> IN (defined) / OUT (changed)
			\@parameter_con,	# List of consensus options											-> IN (defined) / OUT (changed)
			\@parameter_file,	# List of filehandling options										-> IN (defined) / OUT (changed)
			\@parameter_fas,	# List of supermatrix parameter options of FASTA print OUT			-> IN (defined) / OUT (changed)
			\@parameter_3rd,	# List of third position handling									-> IN (defined) / OUT (changed)
			\@parameter_ryc,	# List of RY coding													-> IN (defined) / OUT (changed)
			\@parameter_tra,	# List of Sequence transaltion options								-> IN (defined) / OUT (changed)
			\@parameter_par,	# Print parsimonious sites as extra msa file						-> IN (defined) / OUT (changed)
			\@parameter_ren,	# Rename taxon names of given ifiles								-> IN (defined) / OUT (changed)
			\@parameter_prt,	# Print partition files for concatenated data						-> IN (defined) / OUT (changed)
			\@parameter_pro,	# Start prottest analyses für aa data								-> IN (defined) / OUT (changed)
			\@parameter_mis,	# Replace missing gene sequences by '-' instead using 'X' or 'N'	-> IN (defined) / OUT (changed)
			\%outfile_name,		# Outfilename of output option										-> IN (defined) / OUT (unchanged)
);
##############################





################################################################## END ###############################################################################

sub argv_handling{
	
	my $aref_parameter_all	= $_[0] ;	# List of parameter options of data concatenation					-> IN (defined) / OUT (changed)
	my $aref_parameter_inf	= $_[1] ;	# List of parameter options of info print out						-> IN (defined) / OUT (changed)
	my $aref_parameter_phy	= $_[2] ;	# List of supermatrix parameter options of PHYLIP print OUT			-> IN (defined) / OUT (changed)
	my $aref_parameter_nex	= $_[3] ;	# List of supermatrix parameter options of NEXUS print OUT			-> IN (defined) / OUT (changed)
	my $aref_parameter_con	= $_[4] ;	# List of consensus options											-> IN (defined) / OUT (changed)
	my $aref_parameter_fil	= $_[5] ;	# List of of file handling											-> IN (defined) / OUT (changed)
	my $aref_parameter_fas	= $_[6] ;	# List of supermatrix parameter options of FASTA print OUT			-> IN (defined) / OUT (changed)
	my $aref_parameter_3rd	= $_[7] ;	# List of third position handling									-> IN (defined) / OUT (changed)
	my $aref_parameter_ryc	= $_[8] ;	# List of RY coding													-> IN (defined) / OUT (changed)
	my $aref_parameter_tra	= $_[9] ;	# List of sequence translation options								-> IN (defined) / OUT (changed)
	my $aref_parameter_par	= $_[10];	# List of sequence translation options								-> IN (defined) / OUT (changed)
	my $aref_parameter_ren	= $_[11];	# Rename taxon names of given ifiles								-> IN (defined) / OUT (changed)
	my $aref_parameter_prt	= $_[12];	# Print partition files for concatenated data						-> IN (defined) / OUT (changed)
	my $aref_parameter_pro	= $_[13];	# Start prottest analyses für aa data								-> IN (defined) / OUT (changed)
	my $aref_parameter_mis	= $_[14];	# Replace missing gene sequences by '-' instead using 'X' or 'N'	-> IN (defined) / OUT (changed)
	my $href_outfile_name	= $_[15];	# Outfilename of output option										-> IN (defined) / OUT (unchanged)
	
	
	############################## START ARGV READ IN
	## READ IN and assign terminal command options
	## Set chosen parameter option one step further if associated command option is given
	my ( $commandline ) = join "", @ARGV ;
	
	if ( $commandline ){
		
		
		##############################
		## substitute blanks to '' & split up commandline using '-' signs
		## delete first (empty) entry in @commands
		$commandline =~ s/ |\s+// ;
		my @commands = split '-', $commandline ;
		shift @commands;
		##############################
		
		
		
		##############################
		## Change parameter of given command
		## Skip to main menu if a command option is unknown
		## or '-s' command missing
		for my $single_command ( sort @commands ){
			
			if 		( $single_command =~ /^help$/i	)	{ &help }																	# '-help'	-> open help menu
			elsif 	( $single_command =~ /^i$/i		)	{ @$aref_parameter_inf = reverse @$aref_parameter_inf }						# '-i'		-> set additional info print to NO; 				'-i -i' -> to NO
			elsif 	( $single_command =~ /^f$/i		)	{ @$aref_parameter_all = reverse @$aref_parameter_all }						# '-f'		-> set infile READ IN to defined;					'-f -f' -> to all possible infiles
			elsif 	( $single_command =~ /^a$/i		)	{ @$aref_parameter_fas = reverse @$aref_parameter_fas }						# '-a'		-> set FASTA output to NO;							'-a -a' -> to YES
			elsif 	( $single_command =~ /^d$/i		)	{ @$aref_parameter_3rd = reverse @$aref_parameter_3rd }						# '-d'		-> Reject 3rd position;								'-d -d' -> Remain 3rd position
			elsif 	( $single_command =~ /^j$/i		)	{ @$aref_parameter_par = reverse @$aref_parameter_par }						# '-j'		-> set PARSIMONY output to YES;						'-j -j' -> to NO
			elsif 	( $single_command =~ /^k$/i		)	{ @$aref_parameter_ren = reverse @$aref_parameter_ren }						# '-k'		-> set rename of sequence names to YES;				'-k -k' -> to NO
			elsif 	( $single_command =~ /^l$/i		)	{ @$aref_parameter_prt = reverse @$aref_parameter_prt }						# '-l'		-> set part. file of conc matrix to Supermatrix;	'-l -l' -> to NO
			elsif 	( $single_command =~ /^m$/i		)	{ @$aref_parameter_pro = reverse @$aref_parameter_pro }						# '-m'		-> start prottest under linux;						'-m -m' -> don't start prottest (default)
			elsif 	( $single_command =~ /^g$/i		)	{ @$aref_parameter_mis = reverse @$aref_parameter_mis }						# '-x'		-> replace missing sequences by '-';				'-x -x' -> replace missing sequences by 'X' or 'N' (default)
			elsif 	( $single_command =~ /^e$/i		)	{ my $tl = shift @$aref_parameter_tra ; push @$aref_parameter_tra, $tl }	# '-e'		-> translate nuc data to aa data;					'-e -e' -> translate aa data to nuc data	'-e -e -e' -> no translation
			elsif 	( $single_command =~ /^c$/i		)	{ my $tl = shift @$aref_parameter_con ; push @$aref_parameter_con, $tl }	# '-c'		-> Generate frequency consensus sequence;			'-c -c' -> majority consensus;		 		'-c -c -c' -> strict consensus;	'-c -c -c -c' -> to default (no consensus)
			elsif 	( $single_command =~ /^p$/i		)	{ my $tl = shift @$aref_parameter_phy ; push @$aref_parameter_phy, $tl }	# '-p'		-> set supermatrix format to strict phylip;			'-p -p' -> to relaxed phylip;		 		'-p -p -p' -> to default (no phylip)
			elsif 	( $single_command =~ /^n$/i		)	{ my $tl = shift @$aref_parameter_nex ; push @$aref_parameter_nex, $tl }	# '-n'		-> set supermatrix format to nexus block;			'-n -n' -> to MrBayes nexus block;	 		'-n -n -n' -> to default (no nexus)
			elsif 	( $single_command =~ /^o$/i		)	{ my $tl = shift @$aref_parameter_fil ; push @$aref_parameter_fil, $tl }	# '-o'		-> convert only infiles;							'-o -o' -> Supermatrix & infiles;	 		'-o -o -o' -> print out only supermatrix file
			elsif 	( $single_command =~ /^b$/i		)	{ my $tl = shift @$aref_parameter_ryc ; push @$aref_parameter_ryc, $tl }	# '-b'		-> RY coding complete sequences						'-b -b' -> RY coding 3rd positions;	 		'-b -b -b' -> no RY coding
			elsif 	( $single_command =~ /^s$/i		)	{																			# '-s'		-> Start FASconCAT
				
				&start(
						\@$aref_parameter_all,	# List of parameter options of data concatenation			-> IN (defined) / OUT (changed)
						\@$aref_parameter_inf,	# List of parameter options of info print out				-> IN (defined) / OUT (changed)
						\@$aref_parameter_phy,	# List of supermatrix parameter options of PHYLIP print OUT	-> IN (defined) / OUT (changed)
						\@$aref_parameter_nex,	# List of supermatrix parameter options of NEXUS print OUT	-> IN (defined) / OUT (changed)
						\@$aref_parameter_con,	# List of consensus options									-> IN (defined) / OUT (changed)
						\@$aref_parameter_fil,	# List of of file handling									-> IN (defined) / OUT (changed)
						\@$aref_parameter_fas,	# List of supermatrix parameter options of FASTA print OUT	-> IN (defined) / OUT (changed)
						\@$aref_parameter_3rd,	# List of third position handling							-> IN (defined) / OUT (changed)
						\@$aref_parameter_ryc,	# List of RY coding											-> IN (defined) / OUT (changed)
						\@$aref_parameter_tra,	# List of sequence translation options						-> IN (defined) / OUT (changed)
						\@$aref_parameter_par,	# Print parsimonious sites as extra msa file				-> IN (defined) / OUT (changed)
						\@$aref_parameter_ren,	# Rename taxon names of given ifiles						-> IN (defined) / OUT (changed)
						\@$aref_parameter_prt,	# Print partition files for concatenated data				-> IN (defined) / OUT (changed)
						\@$aref_parameter_pro,	# Start Prottest											-> IN (defined) / OUT (changed)
						\@$aref_parameter_mis,	# Replacement code of missing sequences						-> IN (defined) / OUT (changed)
						\%$href_outfile_name,	# Outfilename of output option								-> IN (defined) / OUT (unchanged)
				)
			}
			else {
				
				print "\n\tCOMMAND-ERROR!: Unknown command ", $single_command, "!\n" ; # 'unknown command'	-> print ERROR prompt
				
				
				
				##############################
				## Open Menu if unknown command is given
				&menu ;
				&parameter (
							
						\@$aref_parameter_all,	# List of parameter options of data concatenation			-> IN (defined) / OUT (changed)
						\@$aref_parameter_inf,	# List of parameter options of info print out				-> IN (defined) / OUT (changed)
						\@$aref_parameter_phy,	# List of supermatrix parameter options of PHYLIP print OUT	-> IN (defined) / OUT (changed)
						\@$aref_parameter_nex,	# List of supermatrix parameter options of NEXUS print OUT	-> IN (defined) / OUT (changed)
						\@$aref_parameter_con,	# List of consensus options									-> IN (defined) / OUT (changed)
						\@$aref_parameter_fil,	# List of of file handling									-> IN (defined) / OUT (changed)
						\@$aref_parameter_fas,	# List of supermatrix parameter options of FASTA print OUT	-> IN (defined) / OUT (changed)
						\@$aref_parameter_3rd,	# List of third position handling							-> IN (defined) / OUT (changed)
						\@$aref_parameter_ryc,	# List of RY coding											-> IN (defined) / OUT (changed)
						\@$aref_parameter_tra,	# List of sequence translation options						-> IN (defined) / OUT (changed)
						\@$aref_parameter_par,	# Print parsimonious sites as extra msa file				-> IN (defined) / OUT (changed)
						\@$aref_parameter_ren,	# Rename taxon names of given ifiles						-> IN (defined) / OUT (changed)
						\@$aref_parameter_prt,	# Print partition files for concatenated data				-> IN (defined) / OUT (changed)
						\@$aref_parameter_pro,	# Start Prottest											-> IN (defined) / OUT (changed)
						\@$aref_parameter_mis,	# Replacement code of missing sequences						-> IN (defined) / OUT (changed)
						\%$href_outfile_name,	# Outfilename of output option								-> IN (defined) / OUT (unchanged)
				)
				##############################
			}
		}
		
		
		
		##############################
		## Open Menu if '-s' command is missing
		&menu ;
		&parameter (
					
					\@$aref_parameter_all,	# List of parameter options of data concatenation			-> IN (defined) / OUT (changed)
					\@$aref_parameter_inf,	# List of parameter options of info print out				-> IN (defined) / OUT (changed)
					\@$aref_parameter_phy,	# List of supermatrix parameter options of PHYLIP print OUT	-> IN (defined) / OUT (changed)
					\@$aref_parameter_nex,	# List of supermatrix parameter options of NEXUS print OUT	-> IN (defined) / OUT (changed)
					\@$aref_parameter_con,	# List of consensus options									-> IN (defined) / OUT (changed)
					\@$aref_parameter_fil,	# List of of file handling									-> IN (defined) / OUT (changed)
					\@$aref_parameter_fas,	# List of supermatrix parameter options of FASTA print OUT	-> IN (defined) / OUT (changed)
					\@$aref_parameter_3rd,	# List of third position handling							-> IN (defined) / OUT (changed)
					\@$aref_parameter_ryc,	# List of RY coding											-> IN (defined) / OUT (changed)
					\@$aref_parameter_tra,	# List of sequence translation options						-> IN (defined) / OUT (changed)
					\@$aref_parameter_par,	# Print parsimonious sites as extra msa file				-> IN (defined) / OUT (changed)
					\@$aref_parameter_ren,	# Rename taxon names of given ifiles						-> IN (defined) / OUT (changed)
					\@$aref_parameter_prt,	# Print partition files for concatenated data				-> IN (defined) / OUT (changed)
					\@$aref_parameter_pro,	# Start Prottest											-> IN (defined) / OUT (changed)
					\@$aref_parameter_mis,	# Replacement code of missing sequences						-> IN (defined) / OUT (changed)
					\%$href_outfile_name,	# Outfilename of output option								-> IN (defined) / OUT (unchanged)
		)
		##############################
	}
	############################## END ARGV READ IN
}

sub menu{ system('cls');
	
	print "";
	printf "\n%68s\n","------------------------------------------------------------"  ;
	printf "%53s\n"  , "Welcome to FASconCAT-G v1.0 !"                                ;
	printf "%58s\n"  , "A perlscript for sequence concatenation"                      ;
	printf "%59s\n"  , "written by Patrick Kueck (ZFMK Bonn, 2010/14)"                ;
	printf "%68s\n\n", "------------------------------------------------------------" ;
}

sub parameter{
	
	my $aref_parameter_all	= $_[0] ;	# List of parameter options of data concatenation			-> IN (defined) / OUT (changed)
	my $aref_parameter_inf	= $_[1] ;	# List of parameter options of info print out				-> IN (defined) / OUT (changed)
	my $aref_parameter_phy	= $_[2] ;	# List of supermatrix parameter options of PHYLIP print OUT	-> IN (defined) / OUT (changed)
	my $aref_parameter_nex	= $_[3] ;	# List of supermatrix parameter options of NEXUS print OUT	-> IN (defined) / OUT (changed)
	my $aref_parameter_con	= $_[4] ;	# List of consensus options									-> IN (defined) / OUT (changed)
	my $aref_parameter_fil	= $_[5] ;	# List of of file handling									-> IN (defined) / OUT (changed)
	my $aref_parameter_fas	= $_[6] ;	# List of supermatrix parameter options of FASTA print OUT	-> IN (defined) / OUT (changed)
	my $aref_parameter_3rd	= $_[7] ;	# List of third position handling							-> IN (defined) / OUT (changed)
	my $aref_parameter_ryc	= $_[8] ;	# List of RY coding											-> IN (defined) / OUT (changed)
	my $aref_parameter_tra	= $_[9] ;	# List of sequence translation options						-> IN (defined) / OUT (changed)
	my $aref_parameter_par	= $_[10];	# List of sequence translation options						-> IN (defined) / OUT (changed)
	my $aref_parameter_ren	= $_[11];	# Rename taxon names of given ifiles						-> IN (defined) / OUT (changed)
	my $aref_parameter_prt	= $_[12];	# Print partition files for concatenated data				-> IN (defined) / OUT (changed)
	my $aref_parameter_pro	= $_[13];	# Start prottest analyses für aa data						-> IN (defined) / OUT (changed)
	my $aref_parameter_mis	= $_[14];	# Replacement code of missing gene sequences 				-> IN (defined) / OUT (changed)
	my $href_outfile_name	= $_[15];	# Outfilename of output option								-> IN (defined) / OUT (unchanged)
	
	
	##############################
	## Print parameter setting and available command options
	print	"\tSTART\t   FASconCAT     :\t\t type <s> <enter>",
			"\n",
			"\n\tINFILES\t   ALL/SINGLE    :\t\t type <f> <enter>",
			"\n\tINFO\t   ALL/BASIC       :\t\t type <i> <enter>",
			"\n",
			"\n\tPROCESSING SMATRIX/CONV. :\t\t type <o> <enter>",
			"\n\tSEQ.TRANS. NO/NUC/AA     :\t\t type <e> <enter>",
			"\n\t3rd POS.   REMAIN/REJECT :\t\t type <d> <enter>",
			"\n\tRY CODING  NO/ALL/3rd    :\t\t type <b> <enter>",
			"\n\tCONSENSUS  N0/YES        :\t\t type <c> <enter>",
			"\n\tRENAME SEQ NO/YES        :\t\t type <k> <enter>",
			"\n\tABSENT SEQ MISSING/INDEL :\t\t type <g> <enter>",
			"\n",
			"\n\tNEXUS      BLOCK/MrBAYES :\t\t type <n> <enter>",
			"\n\tPHYLIP     NO/YES        :\t\t type <p> <enter>",
			"\n\tFASTA      YES/NO        :\t\t type <a> <enter>",
			"\n\tPARTITION  NO/Supermatrix:\t\t type <l> <enter>",
			"\n\tPARSIMONY  NO/YES        :\t\t type <j> <enter>",
			"\n\tPROTTEST   NO/YES        :\t\t type <m> <enter>",
			"\n",
			"\n\tHELP       FASconCAT     :\t\t type <h> <enter>",
			"\n\tQUIT       FASconCAT     :\t\t type <q> <enter>",
			"\n\tPREFACE    FASconCAT     :\t\t type <g> <enter>",
			"\n",
			"\n\t------------------------------------------------------------",
			"\n",
			"\n\tFILE INPUT",
			"\n\t-----------------",
			"\n\tPROCESSING    ALL FILES  :\t$aref_parameter_all->[0]",
			"\n\tPROCESSING SINGLE FILES  :\t$aref_parameter_all->[1]",
			"\n\tPROCESSING SEQ. TRANSL.  :\t$aref_parameter_tra->[0]",
			"\n\tPROCESSING 3rd POSITION\t :\t$aref_parameter_3rd->[0]",
			"\n\tPROCESSING RY CODING\t :\t$aref_parameter_ryc->[0]",
			"\n\tPROCESSING CONSENSUS\t :\t$aref_parameter_con->[0]",
			"\n\tRENAME SEQUENCE NAMES\t :\t$aref_parameter_ren->[0]",
			"\n\tREPLACE ABSENT SEQUENCE\t :\t$aref_parameter_mis->[0]",
			"\n",
			"\n\tFILE OUTPUT",
			"\n\t-----------------",
			"\n\tALL INFO                 :\t$aref_parameter_inf->[0]",
			"\n\tFILE PROCESSING          :\t$aref_parameter_fil->[0]",
			"\n\n\tNEXUS                    :\t$aref_parameter_nex->[0]",
			"\n\tPHYLIP                   :\t$aref_parameter_phy->[0]",
			"\n\tFASTA                    :\t$aref_parameter_fas->[0]\n",
			"\n\tPARTITION FILE           :\t$aref_parameter_prt->[0]",
			"\n\tPARSIMONY                :\t$aref_parameter_par->[0]",
			"\n\tPROTTEST                 :\t$aref_parameter_pro->[0]\n",
			"\n\t------------------------------------------------------------"
	;
	##############################
	
	
	
	##############################
	## (1) Read IN user defined comman
	## (2) If command option allowed, change parameter setting
	## (3) unless user defined command allowed print error prompt
	my			$start_string ;
	my			$single_command	= &commandline ( \$start_string ) ; # (1)
	
	# (2)
	unless (	$single_command	=~ /^s$|^i$|^f$|^q$|^h$|^n$|^p$|^a$|^c$|^o$|^a$|^b$|^q$|^d$|^e$|^g$|^j$|^k$|^l$|^g$|^m$/i ){ print "\n\tCOMMAND-ERROR!: Unknown command ".$single_command."!\n" }
	
	
	# (3)
	if 		( $single_command =~ /^h$/i		)	{ &help }																		# '-h'		-> open help menu
	elsif 	( $single_command =~ /^g$/i		)	{ &preface}																		# '-e'		-> open preface menu
	elsif	( $single_command =~ /^q$/i		)	{ exit }																		# '-q'		-> exit FASconCAT
	elsif 	( $single_command =~ /^i$/i		)	{ @$aref_parameter_inf = reverse @$aref_parameter_inf }							# '-i'		-> set additional info print to YES; 				'-i -i' -> to NO
	elsif 	( $single_command =~ /^f$/i		)	{ @$aref_parameter_all = reverse @$aref_parameter_all }							# '-f'		-> set infile READ IN to defined;					'-f -f' -> to all possible infiles
	elsif 	( $single_command =~ /^a$/i		)	{ @$aref_parameter_fas = reverse @$aref_parameter_fas }							# '-a'		-> set FASTA output to NO;							'-a -a' -> to YES
	elsif 	( $single_command =~ /^j$/i		)	{ @$aref_parameter_par = reverse @$aref_parameter_par }							# '-j'		-> set PARSIMONY output to YES;						'-j -j' -> to NO
	elsif 	( $single_command =~ /^d$/i		)	{ @$aref_parameter_3rd = reverse @$aref_parameter_3rd }							# '-d'		-> Reject 3rd position;								'-d -d' -> Remain 3rd position
	elsif 	( $single_command =~ /^k$/i		)	{ @$aref_parameter_ren = reverse @$aref_parameter_ren }							# '-k'		-> set rename of sequence names to YES;				'-k -k' -> to NO
	elsif 	( $single_command =~ /^l$/i		)	{ @$aref_parameter_prt = reverse @$aref_parameter_prt }							# '-l'		-> set part. file of conc matrix to Supermatrix;	'-l -l' -> to NO
	elsif 	( $single_command =~ /^m$/i		)	{ @$aref_parameter_pro = reverse @$aref_parameter_pro }							# '-m'		-> start prottest under linux;						'-m -m' -> don't start prottest (default)
	elsif 	( $single_command =~ /^g$/i		)	{ @$aref_parameter_mis = reverse @$aref_parameter_mis }							# '-x'		-> replace missing sequences by '-';				'-x -x' -> replace missing sequences by 'X' or 'N' (default)
	elsif 	( $single_command =~ /^e$/i		)	{ my $tl = shift @$aref_parameter_tra ; push @$aref_parameter_tra, $tl }		# '-e'		-> translate nuc data to aa data;					'-e -e' -> translate aa data to nuc data	'-e -e -e' -> no translation
	elsif 	( $single_command =~ /^c$/i		)	{ my $tl = shift @$aref_parameter_con ; push @$aref_parameter_con, $tl }		# '-c'		-> Generate frequency consensus sequence;			'-c -c' -> majority consensus;		 		'-c -c -c' -> strict consensus;	'-c -c -c -c' -> to default (no consensus)
	elsif 	( $single_command =~ /^p$/i		)	{ my $tl = shift @$aref_parameter_phy ; push @$aref_parameter_phy, $tl }		# '-p'		-> set supermatrix format to strict phylip;			'-p -p' -> to relaxed phylip;		 		'-p -p -p' -> to default (no phylip)
	elsif 	( $single_command =~ /^n$/i		)	{ my $tl = shift @$aref_parameter_nex ; push @$aref_parameter_nex, $tl }		# '-n'		-> set supermatrix format to nexus block;			'-n -n' -> to MrBayes nexus block;	 		'-n -n -n' -> to default (no nexus)
	elsif 	( $single_command =~ /^b$/i		)	{ my $tl = shift @$aref_parameter_ryc ; push @$aref_parameter_ryc, $tl }		# '-b'		-> RY coding complete sequences						'-b -b' -> RY coding 3rd positions;	 		'-b -b -b' -> no RY coding
	elsif 	( $single_command =~ /^o$/i		)	{ my $tl = shift @$aref_parameter_fil ; push @$aref_parameter_fil, $tl }		# '-o'		-> convert only infiles;							'-o -o' -> Supermatrix & infiles;	 		'-o -o -o' -> print out only supermatrix file
	elsif 	( $single_command =~ /^s$/i		)	{																				# '-s'		-> Start FASconCAT
		
		&start(
					\@$aref_parameter_all,	# List of parameter options of data concatenation			-> IN (defined) / OUT (changed)
					\@$aref_parameter_inf,	# List of parameter options of info print out				-> IN (defined) / OUT (changed)
					\@$aref_parameter_phy,	# List of supermatrix parameter options of PHYLIP print OUT	-> IN (defined) / OUT (changed)
					\@$aref_parameter_nex,	# List of supermatrix parameter options of NEXUS print OUT	-> IN (defined) / OUT (changed)
					\@$aref_parameter_con,	# List of consensus options									-> IN (defined) / OUT (changed)
					\@$aref_parameter_fil,	# List of of file handling									-> IN (defined) / OUT (changed)
					\@$aref_parameter_fas,	# List of supermatrix parameter options of FASTA print OUT	-> IN (defined) / OUT (changed)
					\@$aref_parameter_3rd,	# List of third position handling							-> IN (defined) / OUT (changed)
					\@$aref_parameter_ryc,	# List of RY coding											-> IN (defined) / OUT (changed)
					\@$aref_parameter_tra,	# List of sequence translation options						-> IN (defined) / OUT (changed)
					\@$aref_parameter_par,	# Print parsimonious sites as extra msa file				-> IN (defined) / OUT (changed)
					\@$aref_parameter_ren,	# Rename taxon names of given ifiles						-> IN (defined) / OUT (changed)
					\@$aref_parameter_prt,	# Print partition files for concatenated data				-> IN (defined) / OUT (changed)
					\@$aref_parameter_pro,	# Start Prottest											-> IN (defined) / OUT (changed)
					\@$aref_parameter_mis,	# Replacement code of missing gene sequences 				-> IN (defined) / OUT (changed)
					\%$href_outfile_name,	# Outfilename of output option								-> IN (defined) / OUT (unchanged)
		)
	}
	
	&menu;
	
	&parameter(
					\@$aref_parameter_all,	# List of parameter options of data concatenation			-> IN (defined) / OUT (changed)
					\@$aref_parameter_inf,	# List of parameter options of info print out				-> IN (defined) / OUT (changed)
					\@$aref_parameter_phy,	# List of supermatrix parameter options of PHYLIP print OUT	-> IN (defined) / OUT (changed)
					\@$aref_parameter_nex,	# List of supermatrix parameter options of NEXUS print OUT	-> IN (defined) / OUT (changed)
					\@$aref_parameter_con,	# List of consensus options									-> IN (defined) / OUT (changed)
					\@$aref_parameter_fil,	# List of of file handling									-> IN (defined) / OUT (changed)
					\@$aref_parameter_fas,	# List of supermatrix parameter options of FASTA print OUT	-> IN (defined) / OUT (changed)
					\@$aref_parameter_3rd,	# List of third position handling							-> IN (defined) / OUT (changed)
					\@$aref_parameter_ryc,	# List of RY coding											-> IN (defined) / OUT (changed)
					\@$aref_parameter_tra,	# List of sequence translation options						-> IN (defined) / OUT (changed)
					\@$aref_parameter_par,	# Print parsimonious sites as extra msa file				-> IN (defined) / OUT (changed)
					\@$aref_parameter_ren,	# Rename taxon names of given ifiles						-> IN (defined) / OUT (changed)
					\@$aref_parameter_prt,	# Print partition files for concatenated data				-> IN (defined) / OUT (changed)
					\@$aref_parameter_pro,	# Start Prottest											-> IN (defined) / OUT (changed)
					\@$aref_parameter_mis,	# Replacement code of missing gene sequences 				-> IN (defined) / OUT (changed)
					\%$href_outfile_name,	# Outfilename of output option								-> IN (defined) / OUT (unchanged)
	) ;
	##############################
}

sub help{ 


	system('cls');

	
print 
<<help
	
	--------------------FASconCAT-G HELP-MENU---------------------------
	
	'Features'
	--------------------------
	- Sequence concatenation of single infiles of different formats:
	  - FASTA   (.fas or .FASTA)
	  - PHYLIP  (.phy -> Relaxed or Strict)
	  - CLUSTAL (.aln)
	  - Individual character fill of absent gene sequences
	    (either by missing chars 'X'|'N' or indels '-')
	- Multiple supermatrix output formats in one process run:
	  - FASTA   (.fas or .FASTA)
	  - PHYLIP  (.phy -> Relaxed or Strict)
	  - NEXUS   (.nex -> BLOCK with or without MrBayes commands)
	- Multiple file conversion of infiles in one process run:
	  - FASTA   (.fas or .FASTA)
	  - PHYLIP  (.phy -> Relaxed or Strict)
	  - NEXUS   (.nex -> BLOCK with or without MrBayes commands)
	- Processing of specific consensus sequences using either...
	  - Strict    consensus rules
	  - Majority  consensus rules or
	  - Frequency consensus rules
	- RY-Coding  of complete sequences  (nucleotide data only)
	- RY-Coding  of 3rd codon positions (nucleotide data only)
	- Removing   of 3rd codon positions (nucleotide data only)
	- Sequence translation (nucleotide data to amino acid data
	  and vice versa)
	- Print OUT of parsimony informative sites within infiles
	  and concatenated supermatrix
	- Handling & identification of specific secondary structure positions
	- Additional output information includes...
	  - Infile ranges in supermatrix sequences
	  - Character state distributions within individual infiles and the supermatrix
	  - Sequence type of individual infiles and the supermatrix
	  - Sequence lengths of individual infiles
	  - Number of gaps and ambiguous sites within infiles and the supermatrix
	  - Number of missing sequences for each individual from the supermatrix 
	  - Number of GC character states within individual infiles and the supermatrix
	  - Number of parsimony informative sites in individual infiles and the supermatrix 
	  - Proportion of character states in supermatrix sequences
	  - Proportion of gaps and ambiguous sites in supermatrix sequences
	  - Proportion of GC in supermatrix sequences
	  - Proportion of parsimony informative sites in supermatrix sequences
	  - Number of structure characters seperated in loops and stems 
	  - Number and percent of loop and stem positions per fragment
	  - Separated list (FcC_structure.txt) for loop positions and 
	    stem pairings whitin the supermatrix and the infiles
	
	- All features can be optionally combined or modified via:
	  - Terminal command line
	  - The FASconCAT (FcC) terminal menu
	--------------------------
	
	
	'Start FASconCAT'
	--------------------------
	To start FASconCAT, open the script by via the terminal.
	
	>Via FASconCAT terminal menu:
	- Open the menu by entering command:
	     perl FASconCAT-G_v1.0.pl <enter> (Linux/Mac)
	          FASconCAT-G_v1.0.pl <enter> (Windows)
	- To start FASconCAT under default, enter command:
	     s <enter>
	
	>Via terminal command line, enter command:
	- perl FASconCAT_v1.0.pl [command options] <enter> (Linux/Mac)
	-      FASconCAT_v1.0.pl [command options] <enter> (Windows)
	- To start FASconCAT under default, enter command:
	     perl FASconCAT_v1.0.pl -s <enter>
	--------------------------
	
	
	Command options:
	--------------------------
	The following command options can be used as stand alone commands
	or in combination with other commands. Command options can be invoked 
	in a single line via the terminal command window or individually via 
	the FASconCAT menu (without the minus'-' specification in the menu). 
	Command options can be input in any order.
	
	
	- Output processing (-o option)
	  - FcC can be used for sequence concatenation and/or file
	    conversion.
	    - Under default, FcC prints OUT a supermatrix with concatenated
	      sequences in FASTA format in a file called
	       "FcC_supermatrix.fas".
	    - For infile format conversions WITHOUT printing
	      a supermatrix, enter command:
	      -o
	    - For format conversions of given infiles AND
	      concatenated supermatrix output, enter command:
	      -o -o
	    
	
	- Sequence translation (-e option)
	  - FcC can translate nucleotide data to amino acid data and 
	    vice versa. Sequence translation happens before exclusion of 
	    third nucleotide codon positions (if -d option is defined).
	  - FcC does not recognize stop codons in amino acid data or 
	    reading frames in nucleotide data. Stop codons have to be 
	    excluded in amino acid data before using FcC!
	  - Nucleotide triplets are translated to their corresponding 
	    amino acid. Ambiguity states in triplets are recognized if 
	    their coding nucleotide states can be definitively assigned 
	    to a specific amino acid state (e.g., 'YTR' -> Leucine/L). 
	    Undefined triplets, such as 'A-G' are translated to '?'. To 
	    select for translation of nucleotide data, enter command:
	    -e
	  - Amino acid states are translated to nucleotide triplets using
	    the compressed IUPAC triplet code for the corresponding amino 
	    acid (e.g. Phenylalanin/F -> 'TTY'). Unrecognized states, such  
	    as '-' or 'X', are translated to '???'.
	    To define reverse translation of amino acid data, enter command:
	    -e -e
	
	- Character fill of lacking gene sequences (-g option)
	  - Under default, lacking sequences in single genes are filled
	    by 'X' (amino-acid) or 'N' (nucleotide): 'MISSING'. With the 
	    -g command, lacking sequences are filled by indel characters
	    '-' instead of missing sequence code: 'INDEL'
	    
	- Renaming sequence names (-k option)
	  - To rename defined sequence names prior to file 
	  	file processing, enter command:
	  	-k
	  -	Note: the user must provide an extra info file named 
	  	"new_seq_names.txt" in the FcC home folder where each row
	  	has the old name delimited from the new name by a tabstop.
	  - Sequences not specified in the "new_seq_names.txt" file are
	  	left unchanged. FcC prints additional information of the 
	  	renaming process to "FcC_rename_control.txt".		    
	
	- Rejecting 3rd codon position (-d option)
	  - To exclude 3rd codon positions in given nucleotide infiles and/or 
	    supermatrices, enter command:
	    -d
	    
	- RY coding (-b option)
	  - To translate third codon positions of nucleotide data into R/Y
	    code (nucleotide data only), enter command:
	    -b
	    
	    e.g., File_1 (nucleotide data):
	    
	    tax1_allel_0 ACGTTTTGGTTT -> ACRTTYTGRTTY...
		
	  - To translate complete sequences into R/Y code, enter command:
	    -b -b
		
	    e.g., File_1 (nucleotide data):
	    
	    tax1_allel_0 ACGTTTTGGTTT -> RYRYYYYRRYYY...     
	
	- Processing consensus sequences (-c option)
	  - FASconCAT can generate consensus sequences for given infile and/or 
	    supermatrix sequences by building consensus states for sequences 
	    with identical sequence names (identical alphanumeric characters 
	    before the first underscore). Consensus outfile sequences are 
	    named by their sequence names, followed by the suffix '_consensus'. 
	    If no underscore is found in sequence names, FcC will use the 
	    complete sequence name for identification of unique sequences.
	    
	    e.g., MSA infile:     ->  Consensus outfile:
	    
	        Taxon_1 AAAACCC
	        Taxon_2 AAAACCC
	        Taxon   AAAACCC ->  Taxon_consensus   AAAACCC
	        
	        taxon1  AAAACCC ->  taxon_1_consensus AAAACCC
	       
	        Tax_1   AAAGCCT
	        Tax_2   AAAGCCT ->  Tax_consensus     AAAGCCT
	
	  - There are three different options for building a consensus sequence:
	    - Most frequent, enter command: -c
	    - Majority rule, enter command: -c -c
	    - Strict, enter command: -c -c -c
	    
	      - Most Frequent ('Freq')
	       - Builds consensus by taking most frequent character state 
	         of each site. If two or more character states are equally frequent,
	         FASconCAT uses the corresponding IUPAC ambiguity code 
	         as the consensus state for nucleotide data or '?' (amino acid data).
			
			 e.g., File_1 (nucleotide data):
			 
			 tax1_allel_0   ACGTTTTCGTTT...
			 tax1_allel_1   ACGTTTTGTTTT...
			 tax1_allel_2   ACGTTTTGGTTT...
			 tax1           ACGTTTTTTTTT...
			 --------------
			 tax1_consensus ACGTTTTGKTTT...
			
			 e.g., File_2 (amino acid data):
			 
			 tax1_allel_0   NYGKRDEDPWFP...
			 tax1_allel_1   NYGKRDEDCWFP...
			 tax1_allel_2   NYGKRDEQCWFP...
			 tax1           NYGKRDEQCWFP...
			 --------------
			 tax1_consensus NYGKRDE?CWFP...
			
	      - Majority Frequent ('Maj'):
	       - Builds consensus by only considering character states that 
	         occur in more than 50% of sequences on a specific site as 
	         consensus state, otherwise assigns '?' (nucleotide & amino 
	         acid data).
			
			 e.g., File_1 (nucleotide or amino acid data):
			 
			 tax1_allel_0   ACGTTTTGGTTT...
			 tax1_allel_1   ACGTTTTGTTTT...
			 tax1_allel_2   ACGTTTTGGTTT...
			 tax1           ACGTTTTTTTTT...
			 --------------
			 tax1_consensus ACGTTTTG?TTT...
			
	      - Strict Consensus ('Strict'):
	       - Builds consensus by taking either sequence states 
	         which are fixed in all sequences at a given site or by 
	         using the corresponding IUPAC ambiguity code (nucleotide
	         data) or 'X' (amino acid data)
			
			 e.g., File_1 (nucleotide data):
			 
			 tax1_allel_0    ACGTTTTGGTTT...
			 tax1_allel_1    ACGTTTTGTTTT...
			 tax1_allel_2    ACGTTTTGGTTT...
			 tax1            ACGTTTTTTTTT...
			 --------------
			 tax1_consensus  ACGTTTTKKTTT...
			 
			 e.g., File_2 (amino acid data):
			 
			 tax1_allel_0   NYGKRDEDPWFP...
			 tax1_allel_1   NYGKRDEDCWFP...
			 tax1_allel_2   NYGKRDEQCWFP...
			 tax1           NYGKRDEQCWFP...
			 --------------
			 tax1_consensus NYGKRDEXXWFP...
		
	- READ IN specific sequence infiles (-f option)
	  - Under default settings, FASconCAT reads in all
	    files in the working directory that are in FASTA
	    (.fas or .FASTA), PHYLIP (.phy) or CLUSTAL (.aln)
	    format. If the '-f' option is chosen, a new list-window 
	    opens after starting FcC. Here input files can be 
	    defined by typing each files assigned numbers (separated
	    by commas and without line spaces) in a single row. 
	    
	    e.g. -f -s <enter> :
	    
	    0 infile_1.fas
	    1 infile_2.phy
	    2 infile_3.aln
	    3 infile_4.FASTA
	    
	    Type: 1,2 <enter> to choose infile_2.phy & infile_3.aln
	    
	- Reduction of FASconCAT sequence information (-i option)
	 - To increase computation speed and reduce the amount of 
	   information FcC prints OUT regarding the infiles and 
	   supermatrix, enter command:
	   -i    
	  
	- Print OUT of parsimony informative sites (-j option)
	 - FASconCAT can print OUT parsimony informative sites
	   for given infiles and the supermatrix (regarded to 
	   the selected -o option). To print OUT parsimony 
	   informative sites, enter command:
	   -j  
	    	
	- Definition of output file formats:
	  - Under default options, FASconCAT prints all output files in 
	    FASTA format. To switch off output of FASTA formatted files, 
	    enter command:
	    -a
	  
	  - FcC can also output NEXUS and PHYLIP formatted files:
	  
	    - To output files in PHYLIP (STRICT) format (sequence 
	      names are restricted to 10 characters), enter command:
	      -p
	      
	    - To output files in PHYLIP (RELAXED) format (sequence names 
	      can be up to 250 characters), enter command:
	      -p -p
	  
	    - To output files in NEXUS (BLOCK) format (NEXUS blocks), 
	      enter command:
	      -n
	      	  
	    - To output files in NEXUS (MrBAYES) format (NEXUS blocks
	      with additional MrBayes commands), enter command:
	      -n -n
	        
	- Print OUT partition file(s) (-l option)
	 - When the sequence concatenation option has been defined, the 
	   user can print additional partition files by entering this
	   command:
	   -l 
	  
	 - If the user selects to output a supermatrix in FASTA and/or 
	   PHYLIP format, FcC will print an associated gene partition file
	   "FcC_supermatrix_partition.txt", which can be directly used for 
	   maximum likelihood analyses in RAxML.
	   
	 - If the user has selected to output a supermatrix in NEXUS (MrBayes) 
	   format, FcC prints an additional NEXUS file 
	   "FcC_supermatrix_partition.nex" in which gene partitions are
	   defined.
	  - Note: the parameters in the partitioned NEXUS (MrBayes) file
	    differ from those in the non-partitioned version.   
	    
	--------------------------
	
	
	For further detailed information please consult 
	the manual or write an email to patrick_kueck\@web.de
	------------------------------------------------------------
	
help
; 

	print  "\tBACK to FASconCAT MAIN-Menu:\t\t type <return>\n"                    ;
	print  "\n\t------------------------------------------------------------\n\t"  ;

	chomp ( my $answer_xy = <STDIN> );

	&menu; &parameter ; 
}

sub preface{

	system('cls');

	
print
<<preface
	
	--------------------FASconCAT PREFACE---------------------
	
	Version     : G-1.0
	Language    : PERL
	Last Update : March, 2021
	Author      : Patrick Kueck, ZFMK Bonn, GERMANY
	e-mail      : patrick_kueck\@web.de
	Homepage    : https://github.com/PatrickKueck/FASconCAT-G
	
	This program is free software; you can distribute it 
	and/or modify it under the terms of the GNU General Public 
	License as published by the Free Software Foundation ; 
	either version 2 of the License, or (at your option) any 
	later version.

	This program is distributed in the hope that it will be 
	useful, but WITHOUT ANY WARRANTY; without even the 
	implied warranty of MERCHANTABILITY or FITNESS FOR A 
	PARTICULAR PURPOSE. See the GNU General Public License for 
	more details. 

	You should have received a copy of the GNU General Public 
	License along with this program; if not, write to the Free 
	Software Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, 
	USA.
	
	For further free downloadable programs visit:
	http://software.zfmk.de
	
	------------------------------------------------------------

preface
; 

	print  "\tBACK to FASconCAT MAIN-Menu:\t\t type <return>\n"                    ;
	print  "\n\t------------------------------------------------------------\n\t"  ;

	chomp ( my $answer_xy = <STDIN> );

	&menu; &parameter ; 
}

sub start{
	
	my $aref_parameter_all	= $_[0] ;	# List of parameter options of data concatenation			-> IN (defined) / OUT (unchanged)
	my $aref_parameter_inf	= $_[1] ;	# List of parameter options of info print out				-> IN (defined) / OUT (unchanged)
	my $aref_parameter_phy	= $_[2] ;	# List of supermatrix parameter options of PHYLIP print OUT	-> IN (defined) / OUT (unchanged)
	my $aref_parameter_nex	= $_[3] ;	# List of supermatrix parameter options of NEXUS print OUT	-> IN (defined) / OUT (unchanged)
	my $aref_parameter_con	= $_[4] ;	# List of consensus options									-> IN (defined) / OUT (unchanged)
	my $aref_parameter_fil	= $_[5] ;	# List of of file handling									-> IN (defined) / OUT (unchanged)
	my $aref_parameter_fas	= $_[6] ;	# List of supermatrix parameter options of FASTA print OUT	-> IN (defined) / OUT (unchanged)
	my $aref_parameter_3rd	= $_[7] ;	# List of third position handling							-> IN (defined) / OUT (unchanged)
	my $aref_parameter_ryc	= $_[8] ;	# List of RY coding											-> IN (defined) / OUT (unchanged)
	my $aref_parameter_tra	= $_[9] ;	# List of sequence translation options						-> IN (defined) / OUT (unchanged)
	my $aref_parameter_par	= $_[10];	# Print parsimonious sites as extra msa file				-> IN (defined) / OUT (unchanged)
	my $aref_parameter_ren	= $_[11];	# Rename taxon names of given ifiles						-> IN (defined) / OUT (unchanged)
	my $aref_parameter_prt	= $_[12];	# Print partition files for concatenated data				-> IN (defined) / OUT (unchanged)
	my $aref_parameter_pro	= $_[13];	# Start prottest analyses for aa data						-> IN (defined) / OUT (unchanged)
	my $aref_parameter_mis	= $_[14];	# Replacement code of missing gene sequences 				-> IN (defined) / OUT (unchanged)
	my $href_outfile_name	= $_[15];	# Outfilename of output option								-> IN (defined) / OUT (unchanged)
	
	
	##############################
	## Print Menu Header and start info
	&menu ;
	
	print "\n\n\t#### FASconCAT-G: START ! ####" ;
	print "\n\t------------------------------------------------------------\n\n" ;
	##############################
	
	
	
	##############################
	## Check if any outputformat has been chosen
	if ( ( $aref_parameter_fas->[0] eq 'NO' ) && ( $aref_parameter_phy->[0] eq 'NO' ) && ( $aref_parameter_nex->[0] eq 'NO' ) ){ die "\n\t!COMMAND-ERROR!: No output format selected!\n" }
	##############################
	
	
	
	##############################
	## Definement of global variabels used in these subroutine and further sub-subroutines
	my (
			%hol_seq_of_tax_of_infile,		# key: input file name; value ; list of taxa and corresponding sequences (e.g. infile_1 = taxon1, sequence_of_taxon1, taxon2, sequence_of_taxon2, taxon3....)
			%taxa_all_files,				# key: taxon name; value : number of occurence in infiles
			%hoh_info_of_infile_of_type,	# key1: infilename, key2: type of info (e.g. 'seqlength'); info type associated value
			@input_files,					# list of input files found in fas, FASTA, phy, or aln format
			$structure_seq					# defined if structure sequence present in infiles
	) ;
	##############################
	
	
	
	##############################
	## READ IN of possible infile names (CLUSTAL, FASTA, PHYLIP formatted files)
	## Store file names in @input_files
	for my $format ( qw/aln phy FASTA fas fasta/ ){
		
		for my $file_input ( <*.$format> ){ push @input_files, $file_input }
	}
	############################## checked !
	
	
	
	##############################
	## If only single, user specified infiles should be concatenated
	## open file request via FASconCAT menu and store only user defined infiles in @input_files
	if ( $aref_parameter_all->[0]  eq 'NO' ){
		
		&single_define (
						\@$aref_parameter_all,	# List of parameter options of data concatenation				-> IN (defined) / OUT (unchanged)
						\@$aref_parameter_inf,	# List of parameter options of info print out					-> IN (defined) / OUT (unchanged)
						\@$aref_parameter_phy,	# List of supermatrix parameter options of PHYLIP print OUT		-> IN (defined) / OUT (unchanged)
						\@$aref_parameter_nex,	# List of supermatrix parameter options of NEXUS print OUT		-> IN (defined) / OUT (unchanged)
						\@$aref_parameter_con,	# List of consensus options										-> IN (defined) / OUT (unchanged)
						\@$aref_parameter_fil,	# List of of file handling										-> IN (defined) / OUT (unchanged)
						\@$aref_parameter_fas,	# List of supermatrix parameter options of FASTA print OUT		-> IN (defined) / OUT (unchanged)
						\@$aref_parameter_3rd,	# List of third position handling								-> IN (defined) / OUT (unchanged)
						\@$aref_parameter_ryc,	# List of RY coding												-> IN (defined) / OUT (unchanged)
						\@$aref_parameter_tra,	# List of sequence translation options							-> IN (defined) / OUT (unchanged)
						\@$aref_parameter_par,	# Print parsimonious sites as extra msa file					-> IN (defined) / OUT (unchanged)
						\@$aref_parameter_ren,	# Rename taxon names of given ifiles							-> IN (defined) / OUT (unchanged)
						\@$aref_parameter_prt,	# Print partition files for concatenated data					-> IN (defined) / OUT (unchanged)
						\@$aref_parameter_pro,	# Start prottest analyses for aa data							-> IN (defined) / OUT (unchanged)
						\@$aref_parameter_mis,	# Replacement code of missing gene sequences 					-> IN (defined) / OUT (unchanged)
						\%$href_outfile_name,	# Outfilename of output option									-> IN (defined) / OUT (unchanged)
						\@input_files			# list of input files found in fas, FASTA, phy, or aln format	-> IN (defined) / OUT (unchanged)
		);
	} 
	##############################
	
	
	
	############################################################ START FILE READ IN & FILE CHECK
	## READ IN file content and check for correct format of infiles				-> &input_check
	&input_check (
					\%hol_seq_of_tax_of_infile,		# key: input file name; value ; list of taxa and corresponding sequences				-> IN (undefined) / OUT (defined)
					\%taxa_all_files,				# key: taxon name; value : number of occurence in infiles								-> IN (undefined) / OUT (defined)
					\@input_files,					# list of input files found in fas, FASTA, phy, or aln format							-> IN (defined) / OUT (unchanged)
					\%hoh_info_of_infile_of_type,	# key1: infilename, key2: type of info (e.g. 'seqlength'); info type associated value	-> IN (undefined) / OUT (defined)
					\$structure_seq				# defined if structure sequence present in infiles											-> IN (undefined) / OUT (defined)
	) ;
	############################################################ END FILE READ IN & FILE CHECK
	
	
	
	############################################################ START SEQ RENAMING
	unless ( $aref_parameter_ren->[0] eq 'NO' ){
		
		&seq_renaming(
					\%hol_seq_of_tax_of_infile,		# key: input file name; value ; list of taxa and corresponding sequences				-> IN (defined) / OUT (changed)
					\%hoh_info_of_infile_of_type,	# key1: infilename, key2: type of info (e.g. 'seqlength'); info type associated value	-> IN (defined) / OUT (changed)
					\%taxa_all_files,				# key: taxon name; value : number of occurence in infiles								-> IN (defined) / OUT (changed)
					\$structure_seq				# defined if structure sequence present in infiles											-> IN (defined) / OUT (changed)
		);
	}
	############################################################ END SEQ RENAMING
	
	
	
	############################################################ START SEQ TRANSLATION
	unless ( $aref_parameter_tra->[0] eq 'NO' ){
		
		&seq_translation(
					\%hol_seq_of_tax_of_infile,		# key: input file name; value ; list of taxa and corresponding sequences				-> IN (defined) / OUT (changed)
					\%hoh_info_of_infile_of_type,	# key1: infilename, key2: type of info (e.g. 'seqlength'); info type associated value	-> IN (defined) / OUT (changed)
					\$aref_parameter_tra->[0]		# defined ry coding parameter option													-> IN (defined) / OUT (unchanged)
		);
	}
	############################################################ END SEQ TRANSLATION
	
	
	
	############################################################ START CONSENSUS PROCESS OF SEQUENCE BLOCKS IF DEFINED
	## If consensus setup is set to yes ('Freq', 'Maj', or 'Strict') built consensus sequence of defined sequence blocks (sequences with identic taxon names before the first underscore)
	## regarded to chosen consensus option and store consensus sequence of defined sequence blocks in...
	## ...%consensus_seq_of_taxon -> key consensus taxon name (original taxon prefix before first underscore with additional '_consensus' suffix); value: consensus sequence
	## -----------------
	## e.g.:
	## tax1_allel_0								->	ACGTTTTCGTTT...
	## tax1_allel_1								->	ACGTTTTAATTT...
	## tax1_allel_2								->	ACGTTTTGCTTT...
	## tax1										->	ACGTTTTTTTTT...
	## $consensus_seq_of_taxon{tax1_consensus}	= 	ACGTTTTNNTTT...
	## -----------------
	## Consensus sequence built in subroutine -> &make_consensus
	unless ( $aref_parameter_con->[0] eq 'NO' ){
		
		
		##############################
		## Create consensus sequences of defined sequence blocks in given infile
		&make_consensus (
							\%hol_seq_of_tax_of_infile,		# key: input file name; value ; list of taxa and corresponding sequences				-> IN (defined) / OUT (changed)
							\$aref_parameter_con->[0],		# defined consensus parameter (Freq, Maj, or Strict)									-> IN (defined) / OUT (unchanged)
							\%hoh_info_of_infile_of_type,	# key1: infilename, key2: type of info (e.g. 'seqlength'); info type associated value	-> IN (defined) / OUT (unchanged)
							\%taxa_all_files,				# key: taxon name; value : number of occurence in infiles								-> IN (defined) / OUT (changed)
		)
		##############################
	}
	############################################################ END CONSENSUS PROCESS OF SEQUENCE BLOCKS IF DEFINED
	
	
	
	############################################################ START RY CODING
	## If defined, RY coding of either complete nucleotide sequences or 3rd positions of nuc sequences
	unless ( $aref_parameter_ryc->[0] eq 'NO' ){
		
		&ry_coding(
					\%hol_seq_of_tax_of_infile,		# key: input file name; value ; list of taxa and corresponding sequences				-> IN (defined) / OUT (changed)
					\%hoh_info_of_infile_of_type,	# key1: infilename, key2: type of info (e.g. 'seqlength'); info type associated value	-> IN (defined) / OUT (unchanged)
					\$aref_parameter_ryc->[0]		# defined ry coding parameter option													-> IN (defined) / OUT (unchanged)
		) ;
	}
	############################################################ END RY CODING
	
	
	
	############################################################ START REMOVING 3RD SEQUENCE POSITION
	## If defined, remove third sequence position (only in nucleotide data)
	unless ( $aref_parameter_3rd->[0] eq 'Remain' ){
		
		&reject_third_nuc_position(
										\%hol_seq_of_tax_of_infile,		# key: input file name; value ; list of taxa and corresponding sequences				-> IN (defined) / OUT (changed)
										\%hoh_info_of_infile_of_type,	# key1: infilename, key2: type of info (e.g. 'seqlength'); info type associated value	-> IN (defined) / OUT (changed)
		) ;
	}
	############################################################ START REMOVING 3RD SEQUENCE POSITION
	
	
	
	########################################################### START SEQUENCE CONCATENATION
	## Sequence concatenation and extraction of basal sequence information
	## unless file handling parameter is set to 'Convert' -> Convert single files to other output format without adding missing sequences
	unless ( $aref_parameter_fil->[0]  eq 'Convert' ){
		
		&concatenate (
						
						\%hol_seq_of_tax_of_infile,		# key: input file name; value ; list of taxa and corresponding sequences				-> IN (defined) / OUT (changed)
						\%hoh_info_of_infile_of_type,	# key1: infilename, key2: type of info (e.g. 'seqlength'); info type associated value	-> IN (undefined) / OUT (defined)
						\%taxa_all_files,				# key: taxon name; value : number of occurence in infiles								-> IN (defined) / OUT (unchanged)
						\$structure_seq,				# defined if structure sequence present in infiles										-> IN (defined) / OUT (unchanged)
						\$aref_parameter_mis->[0],		# defined replacement parameter (MISSING or INDEL)										-> IN (defined) / OUT (unchanged)
		)
	}
	########################################################### END SEQUENCE CONCATENATION
	
	
	
	############################################################ START PROTTEST
	## If defined, RY coding of either complete nucleotide sequences or 3rd positions of nuc sequences
	unless ( $aref_parameter_pro->[0] eq 'NO' ){
		
		unless ( ( $aref_parameter_fas->[0] eq 'NO' ) && ( $aref_parameter_phy->[0] eq 'NO' ) ){
			
			&start_prottest(
						
						\%hol_seq_of_tax_of_infile,			# key: input file name; value ; list of taxa and corresponding sequences				-> IN (defined) / OUT (changed)
						\%hoh_info_of_infile_of_type,		# key1: infilename, key2: type of info (e.g. 'seqlength'); info type associated value	-> IN (defined) / OUT (unchanged)
						\$aref_parameter_fil->[0],			# List of of file handling																-> IN (defined) / OUT (unchanged)
						\$aref_parameter_prt->[0],			# Print partition files for concatenated data											-> IN (defined) / OUT (unchanged)
						\$href_outfile_name->{prottest},	# PROTTEST SCRIPT NAME																	-> IN (defined) / OUT (unchanged)
			) ;
		}
		else{ print "\n\n\t!COMMAND-ERROR!: ProtTest not started!\n\tFasta or phylip outfile format must be defined for ProtTest analyses!\n" }
	}
	############################################################ END PROTTEST
	
	
	
	##############################
	## If defined nexus output set additional info option to 'all info' -> 'YES'
	if ( $aref_parameter_nex->[0]  =~ /BLOCK|MrBAYES/ ){ $aref_parameter_inf->[0] = 'YES' }
	##############################
	
	
	
	############################################################ START STRUCTURE SEQUENCE ANALYSIS
	## If additional parameter info option set to yes and defined structure sequence found in infiles
	## extract structure sequence info like loop and stem positions; N of states...
	if ( ( $aref_parameter_inf->[0] eq 'YES' ) && ( $structure_seq ) ){
		
		&structure_handling (
						
						\%hol_seq_of_tax_of_infile,		# key: input file name; value ; list of taxa and corresponding sequences				-> IN (defined) / OUT (unchanged)
						\%hoh_info_of_infile_of_type,	# key1: infilename, key2: type of info (e.g. 'seqlength'); info type associated value	-> IN (defined) / OUT (changed)
		) ;
	}
	############################################################ END STRUCTURE SEQUENCE ANALYSIS
	
	
	
	############################################################ START EXTRACTION OF ADDITIONAL SEQUENCE STATE INFO
	## Count additional info about number of character states for each file sorted by missing data, ambiguity, indel event, and informative state
	## if additional info options is set to 'YES'
	if ( $aref_parameter_inf->[0] eq "YES" ){
		
		&get_info (
						
						\%hol_seq_of_tax_of_infile,			# key: input file name; value ; list of taxa and corresponding sequences				-> IN (defined) / OUT (unchanged)
						\%hoh_info_of_infile_of_type,		# key1: infilename, key2: type of info (e.g. 'seqlength'); info type associated value	-> IN (defined) / OUT (changed)
		) ;
	}
	############################################################ END EXTRACTION OF ADDITIONAL SEQUENCE STATE INFO
	
	
	
	############################################################ START PRINT OUT OF SUPERMATRIX & EXTRACTED INFO
	## Print OUT outfiles
	&print_out (
						
						\%hol_seq_of_tax_of_infile,			# key: input file name; value ; list of taxa and corresponding sequences				-> IN (defined) / OUT (unchanged)
						\%hoh_info_of_infile_of_type,		# key1: infilename, key2: type of info (e.g. 'seqlength'); info type associated value	-> IN (defined) / OUT (unchanged)
						\%$href_outfile_name,				# Outfilename of output option															-> IN (defined) / OUT (unchanged)
						\@$aref_parameter_all,				# List of parameter options of data concatenation										-> IN (defined) / OUT (unchanged)
						\@$aref_parameter_inf,				# List of parameter options of info print out											-> IN (defined) / OUT (unchanged)
						\@$aref_parameter_phy,				# List of parameter options of PHYLIP print OUT											-> IN (defined) / OUT (unchanged)
						\@$aref_parameter_nex,				# List of supermatrix parameter options of NEXUS print OUT								-> IN (defined) / OUT (unchanged)
						\@$aref_parameter_con,				# List of consensus options																-> IN (defined) / OUT (unchanged)
						\@$aref_parameter_fil,				# List of filehandling options															-> IN (defined) / OUT (unchanged)
						\@$aref_parameter_fas,				# List of parameter options of FASTA print OUT											-> IN (defined) / OUT (unchanged)
						\@$aref_parameter_3rd,				# List of third position handling														-> IN (defined) / OUT (unchanged)
						\@$aref_parameter_ryc,				# List of RY coding																		-> IN (defined) / OUT (unchanged)
						\@$aref_parameter_tra,				# List of sequence translation options													-> IN (defined) / OUT (unchanged)
						\@$aref_parameter_par,				# Print parsimonious sites as extra msa file											-> IN (defined) / OUT (unchanged)
						\@$aref_parameter_ren,				# Rename taxon names of given ifiles													-> IN (defined) / OUT (unchanged)
						\@$aref_parameter_prt,				# Print partition files for concatenated data											-> IN (defined) / OUT (unchanged)
						\@$aref_parameter_pro,				# Start prottest analyses for aa data													-> IN (defined) / OUT (unchanged)
						\@$aref_parameter_mis,				# Replacement code of missing gene sequences 											-> IN (defined) / OUT (unchanged)
						\$structure_seq						# defined if structure sequence present in infiles										-> IN (defined) / OUT (unchanged)
	
	) ;
	############################################################ END PRINT OUT OF SUPERMATRIX & EXTRACTED INFO
	
	
	
	############################################################ START PRINT OUT OF SUPERMATRIX & EXTRACTED INFO
	## Print END section on terminal
	&end (
						
						\%$href_outfile_name,				# Outfilename of output option															-> IN (defined) / OUT (unchanged)
						\@$aref_parameter_inf,				# List of parameter options of info print out											-> IN (defined) / OUT (unchanged)
						\@$aref_parameter_phy,				# List of parameter options of PHYLIP print OUT											-> IN (defined) / OUT (unchanged)
						\@$aref_parameter_nex,				# List of supermatrix parameter options of NEXUS print OUT								-> IN (defined) / OUT (unchanged)
						\@$aref_parameter_con,				# List of consensus options																-> IN (defined) / OUT (unchanged)
						\@$aref_parameter_fil,				# List of filehandling options								-							-> IN (defined) / OUT (unchanged)
						\@$aref_parameter_fas,				# List of parameter options of FASTA print OUT											-> IN (defined) / OUT (unchanged)
						\@$aref_parameter_3rd,				# List of third position handling														-> IN (defined) / OUT (unchanged)
						\@$aref_parameter_ryc,				# List of RY coding																		-> IN (defined) / OUT (unchanged)
						\@$aref_parameter_tra,				# List of Sequence transaltion options													-> IN (defined) / OUT (unchanged)
						\@$aref_parameter_par,				# Print parsimonious sites as extra msa file											-> IN (defined) / OUT (unchanged)
						\@$aref_parameter_ren,				# Rename taxon names of given ifiles													-> IN (defined) / OUT (unchanged)
						\@$aref_parameter_prt,				# Print partition files for concatenated data											-> IN (defined) / OUT (unchanged)
						\@$aref_parameter_mis,				# Replacement code of missing gene sequences 											-> IN (defined) / OUT (unchanged)
						\$structure_seq						# defined if structure sequence present in infiles										-> IN (defined) / OUT (unchanged)
	) ;
	############################################################ END PRINT OUT OF SUPERMATRIX & EXTRACTED INFO
	
	sub start_prottest{
		
		my $href_hol_seq_of_tax_of_infile		= $_[0] ;	# key: input file name; value ; list of taxa and corresponding sequences				-> IN (defined) / OUT (unchanged)
		my $href_hoh_info_of_infile_of_type		= $_[1] ;	# key1: infilename, key2: type of info (e.g. 'seqlength'); info type associated value	-> IN (defined) / OUT (unchanged)
		my $sref_file_parameter					= $_[2] ;	# file handling option																	-> IN (defined) / OUT (unchanged)
		my $sref_part_parameter					= $_[3] ;	# partition file handling																-> IN (defined) / OUT (unchanged)
		my $sref_prottest_name					= $_[4] ;	# Outfilename of output option															-> IN (defined) / OUT (unchanged)
		
		print "\n\tStart PhyML-ProtTest..."; 
		
		
		
		##############################
		## ProtTest will only be started if sequence concatenation and partition file print out have been selected
		unless	( $$sref_file_parameter =~ /Supermatrix/ ){ print "\n\t!COMMAND-ERROR!: ProtTest not started due to unspecified concatenation process!" }
		elsif 	( $$sref_part_parameter =~ /NO/ ){ print "\n\t!COMMAND-ERROR!: ProtTest not started due to unspecified partition file print out !" }
		else{
			
			######################################################## START PROTTEST RUN FOR SINGLE INFILES
			PROTFILE:
			for my $infile ( sort {$a<=>$b} keys %$href_hol_seq_of_tax_of_infile ){
				
				
				
				############################################################ START HANDLING OF NUC DATA
				# PROTTEST cannot be used for nucleotide data....
				if ( $href_hoh_info_of_infile_of_type->{$infile}{seqtype} eq 'nu' ){
					
					
					
					##############################
					## Assign DNA model to given infile for RAxML partition file
					$href_hoh_info_of_infile_of_type->{$infile}{protmodel} = 'DNA' ;
					##############################
					
					
					
					##############################
					## Print file handling info
					print "\n\tCannot Start PROTTEST of NUC file ", $infile, "..." ;
					##############################
				}
				############################################################ END HANDLING OF NUC DATA
				
				
				
				############################################################ START HANDLING OF AA DATA
				else{
					
					
					##############################
					## Prottest must not be conducted for concatenated supermatrix
					if ( $infile eq 'supermatrix' ){ next PROTFILE }
					##############################
					
					
					
					##############################
					## print file handling info
					( my $file = $infile ) =~ s/.fas$|.fasta$|.FASTA$|.phy$|.aln$// ;
					print "\n\t...for ", $file, "\n\n" ; 
					##############################
					
					
					
					##############################
					## PRINT out temporary phylip format for prottest
					my @data	= exists( $href_hol_seq_of_tax_of_infile->{$infile}) ? @{$href_hol_seq_of_tax_of_infile->{$infile}} :( ) ;
					
					my $Ntaxa	= @data / 2 ;
					my $Nchar	= length $data[1] ;
					
					open  OUTphy,	">prottest_infile.tmp" || warn "\n\t!FILE-ERROR!: Cannot OPEN OUT 'prottest_infile.tmp'!\n" ;
					print OUTphy	"$Ntaxa $Nchar\n" ;
					
					while ( @data ){
							
							my	$taxon				= shift @data ;
								$taxon				=~ s/ /_/g;
							my	$sequence_of_taxon	= shift @data ;
							
							print OUTphy $taxon, "   ", $sequence_of_taxon, "\n"
					}
					##############################
					
					
					
					##############################
					## START PROTTEST
					( my $file = $infile ) =~ s/.fas$|.FAS$|.FASTA$|.phy$|.aln$|.fasta$// ;
					my $prottest_out = "FcC_prottest_results_".$file.".txt" ;
					system ( "java -jar $$sref_prottest_name -i prottest_infile.tmp -o $prottest_out" ) ;
					##############################
					
					
					
					##############################
					## Delete prottest temp infile
					unlink ( 'prottest_infile.tmp' );
					##############################
					
					
					
					##############################
					## Read IN of the prottest result file and extraction of model with best BIC criterium
					open INpt, "<$prottest_out" || warn "\n\t!FILE-ERROR!: Cannot READ IN ", $prottest_out , "'!\n" ;
					while ( my $line = <INpt> ){ chomp $line;
						
						if ( $line =~ /^Best model according to BIC: \w+$/ ){
							
							my @prts			= split ": ", $line ;
							my $best_bic_model	= $prts[1] ;
							#print "\n", $best_bic_model, "\n" ;
							
							
							##############################
							## Assign best BIC model to given infile for RAxML partition file
							$href_hoh_info_of_infile_of_type->{$infile}{protmodel} = $best_bic_model ;
							##############################
						}
					}
					##############################
				}
			}
		}
		############################################################ END PROTTEST RUN FOR SINGLE INFILES
	}
	
	sub seq_renaming{
		
		my $href_hol_seq_of_tax_of_infile		= $_[0] ;	# key: input file name; value ; list of taxa and corresponding sequences				-> IN (defined) / OUT (changed)
		my $href_hoh_info_of_infile_of_type		= $_[1] ;	# key1: infilename, key2: type of info (e.g. 'seqlength'); info type associated value	-> IN (defined) / OUT (changed)
		my $href_name_list						= $_[2] ;	# key: taxon name; value : number of occurence in infiles								-> IN (defined) / OUT (changed)
		my $sref_structure						= $_[3] ;	# defined if structure sequence present in infiles										-> IN (defined) / OUT (changed)
		
		
		########################################################################################## START RENAMING OF INFILES
		
		
		
		##############################
		## READ IN file with new taxon names
		my $codefile = "new_seq_names.txt" ;
		open INcode, "<$codefile" || die "\n\t!FILE-ERROR! Cannot READ IN codefile ", $codefile, "!\n" ;
		
		print "\n\tREAD IN codefile ", $codefile, " with new sequence names..." ;
		####################################
		
		
		
		####################################
		## OPEN OUT Control.txt file
		open OUTtxt, ">FcC_rename_control.txt" || die "\n\t!FILE-ERROR! Cannot OPEN OUT control file FcC_rename_control.txt!\n" ;
		####################################
		
		
		my (
				%new_name_of_old_name,	# key: old sequence name; value: new sequence name
				%defined_new_name
		);
		
		
		####################################
		## read in each line of defined clade file
		## store new sequence name (second row) as hash value of old sequence name (first row)
		## if no new sequence name dafinde take old sequence name as new sequence name
		## count number of defined new names for each old sequence name
		my $linecounter = 1 ;
		while ( my $line = <INcode> ){
			
			
			####
			# exclude newline signs of given line and break line by '\t'
			# die with an error prompt if first element of @parts (old names) is empty
			chomp $line; #$line =~ s/\s+//g;
			my @parts = split "\t", $line ;
			unless ( @parts == 2			){ die "\n\t!FILE-ERROR! more as 1 tabstop sign found in new_seq_names.txt\n\tBe aware of correct file format\n" }
			#unless ( $parts[0] =~ /^\w+$/ 	){ die "\n\t!FILE-ERROR! Unallowed sign(s) in old sequence name in line ", $linecounter, " of codefile ", $codefile, "!\n\tBe aware of correct file format\n" }
			#unless ( $parts[1] =~ /^\w+$/ 	){ die "\n\t!FILE-ERROR! Unallowed sign(s) in new sequence name in line ", $linecounter, " of codefile ", $codefile, "!\n\tBe aware of correct file format\n" }
			####
			
			unless ( $defined_new_name{$parts[0]} ){ $defined_new_name{$parts[0]}++ } else{ die "\n\t!FILE-ERROR! New sequence name ", $parts[0], " defined multiple times in ", $codefile, "!\n" }
			
			####
			# if new name dfined store new name as hash value of old name
			# count total number of defined new names for old name
			$new_name_of_old_name{$parts[0]} = $parts[1]  ;
			print OUTtxt	"\n\tFound new name\tfor ", $parts[0], "\t->\t", $parts[1] ;
			####
			
			
			
			####
			# set line counter +1
			$linecounter++
			####
		}
		
		print OUTtxt	"\n";
		close INcode ;
		####################################
		
		
		
		####################################
		## empty old sequence name list of all infiles
		%$href_name_list = () ;
		####################################
		
		
		
		####################################
		## READ IN FASTA file
		## Define new fasta outfile ($fasta_new)
		## Change sequence names, Print OUT sequences and changed sequence names in FASTA format
		for my $infile ( sort {$a<=>$b} keys %$href_hol_seq_of_tax_of_infile ){
			
			
			##############################
			## print info
			print 			"\n\tRenaming sequences of ", $infile, "..." ;
			print OUTtxt 	"\n\tRenaming sequences of ", $infile, "..." ;
			##############################
			
			
			
			##############################
			## Store identified sequences and associated taxon names of given infile in...
			## ...%sequence_of_taxon -> key: original taxon name; value: associated sequence
			my $aref_list_taxon_sequence_of_file	= exists( $href_hol_seq_of_tax_of_infile->{$infile}) ? \@{$href_hol_seq_of_tax_of_infile->{$infile}} :( ) ;
			my %sequence_of_taxon						= @$aref_list_taxon_sequence_of_file ;
			##############################
			
			
			##############################
			## START RENAMING OF SEQUENCES OF GIVEN INFILE
			my $counter_changed		= 0 ;
			my $counter_unchanged	= 0 ;
			my %sequence_of_new_names ;
			for my $taxon ( sort keys %sequence_of_taxon ){
				
				
				if ( $new_name_of_old_name{$taxon} ){
					
					$sequence_of_new_names{$new_name_of_old_name{$taxon}} = $sequence_of_taxon{$taxon} ;
					$counter_changed++ ;
					
					print OUTtxt	"\n\t", $taxon, "\tchanged to\t", $new_name_of_old_name{$taxon} ;
					
					
					
					##############################
					## change structure name of given variabels
					if ( $taxon eq $$sref_structure ){
						
						$$sref_structure = $new_name_of_old_name{$taxon} ;
						$href_hoh_info_of_infile_of_type->{$infile}{taxstruct} = $$sref_structure
					}
					##############################
				}
				
				else{
					
					$sequence_of_new_names{$taxon} = $sequence_of_taxon{$taxon} ;
					$counter_unchanged++;
					
					print OUTtxt	"\n\t", $taxon, "\tnot changed"
				}
			}
			##############################
			
			
			
			##############################
			## store new sequence nnames for given infiles and in the overall haslist of sequence names
			$href_hol_seq_of_tax_of_infile->{$infile} = () ;
			
			for my $taxon ( keys %sequence_of_new_names ){
				
				push @{$href_hol_seq_of_tax_of_infile->{$infile}}, ( $taxon, $sequence_of_new_names{$taxon} ) ;
				$href_name_list->{$taxon}++
			}
			
			%sequence_of_new_names = () ;
			##############################
			
			
			
			##############################
			## info print
			print OUTtxt "\n\t", $infile, ":\t", $counter_changed, " sequence names changed\t", $counter_unchanged, " sequence names not changed\n" ;
			##############################
		}
		
		close OUTtxt;
		########################################################################################## END RENAMING OF INFILES
	}
	
	sub single_define{
		
		my $aref_parameter_all	= $_[0] ;	# List of parameter options of data concatenation				-> IN (defined) / OUT (unchanged)
		my $aref_parameter_inf	= $_[1] ;	# List of parameter options of info print out					-> IN (defined) / OUT (unchanged)
		my $aref_parameter_phy	= $_[2] ;	# List of supermatrix parameter options of PHYLIP print OUT		-> IN (defined) / OUT (unchanged)
		my $aref_parameter_nex	= $_[3] ;	# List of supermatrix parameter options of NEXUS print OUT		-> IN (defined) / OUT (unchanged)
		my $aref_parameter_con	= $_[4] ;	# List of consensus options										-> IN (defined) / OUT (unchanged)
		my $aref_parameter_fil	= $_[5] ;	# List of of file handling										-> IN (defined) / OUT (unchanged)
		my $aref_parameter_fas	= $_[6] ;	# List of supermatrix parameter options of FASTA print OUT		-> IN (defined) / OUT (unchanged)
		my $aref_parameter_3rd	= $_[7] ;	# List of third position handling								-> IN (defined) / OUT (unchanged)
		my $aref_parameter_ryc	= $_[8] ;	# List of RY coding												-> IN (defined) / OUT (unchanged)
		my $aref_parameter_tra	= $_[9] ;	# List of sequence translation options							-> IN (defined) / OUT (unchanged)
		my $aref_parameter_par	= $_[10];	# Print parsimonious sites as extra msa file					-> IN (defined) / OUT (unchanged)
		my $aref_parameter_ren	= $_[11];	# Rename taxon names of given ifiles							-> IN (defined) / OUT (unchanged)
		my $aref_parameter_prt	= $_[12];	# Print partition files for concatenated data					-> IN (defined) / OUT (unchanged)
		my $aref_parameter_pro	= $_[13];	# Start prottest analyses for aa data							-> IN (defined) / OUT (unchanged)
		my $aref_parameter_mis	= $_[14];	# Replacement code of missing gene sequences 					-> IN (defined) / OUT (unchanged)
		my $href_outfile_name	= $_[15];	# Outfilename of output option									-> IN (defined) / OUT (unchanged)
		my $aref_input_files	= $_[16];	# list of input files found in fas, FASTA, phy, or aln format	-> IN (defined) / OUT (changed)
		
		
		##############################
		## Print possible infiles (sorted by name) with associated filenumber (array number)
		## -----------------
		## e.g.:
		## 0	infile1.fas
		## 1	infile2.phy
		## 2	infile3.aln
		## -----------------
		my $infilelist_counter = 0 ; @$aref_input_files = sort @$aref_input_files ;
		for ( @$aref_input_files ){ print "\n\t", $infilelist_counter, "\t", $_ ; $infilelist_counter++ }
		
		$infilelist_counter = () ;
		##############################
		
		
		
		##############################
		## Print info text line, start command request -> &commandline
		## If user defined command equal 'b' || 'B' -> go back to main menu -> &menu; &parameter
		## else check for correct command format and repeat command request until correct format is given
		## -----------------
		## correct chosen input file examples:
		## 0,2
		## 1
		## 0,1,2
		## -----------------
		my		$info_text		=  "\n\tNumber of INPUT files comma separated ( BACK <b> ): " ;
		my		$number_input	=  () ;
		until ( $number_input	=~ /^\d+(,\d+)*$|^b$|^B$/i ){ $number_input = &commandline( \$info_text ) }
				$number_input	=~ /^b$/i and do { 
					
					&menu; 
					&parameter(
								\@$aref_parameter_all,	# List of parameter options of data concatenation				-> IN (defined) / OUT (unchanged)
								\@$aref_parameter_inf,	# List of parameter options of info print out					-> IN (defined) / OUT (unchanged)
								\@$aref_parameter_phy,	# List of supermatrix parameter options of PHYLIP print OUT		-> IN (defined) / OUT (unchanged)
								\@$aref_parameter_nex,	# List of supermatrix parameter options of NEXUS print OUT		-> IN (defined) / OUT (unchanged)
								\@$aref_parameter_con,	# List of consensus options										-> IN (defined) / OUT (unchanged)
								\@$aref_parameter_fil,	# List of of file handling										-> IN (defined) / OUT (unchanged)
								\@$aref_parameter_fas,	# List of supermatrix parameter options of FASTA print OUT		-> IN (defined) / OUT (unchanged)
								\@$aref_parameter_3rd,	# List of third position handling								-> IN (defined) / OUT (unchanged)
								\@$aref_parameter_ryc,	# List of RY coding												-> IN (defined) / OUT (unchanged)
								\@$aref_parameter_tra,	# List of sequence translation options							-> IN (defined) / OUT (unchanged)
								\@$aref_parameter_par,	# Print parsimonious sites as extra msa file					-> IN (defined) / OUT (unchanged)
								\@$aref_parameter_ren,	# Rename taxon names of given ifiles							-> IN (defined) / OUT (unchanged)
								\@$aref_parameter_prt,	# Print partition files for concatenated data					-> IN (defined) / OUT (unchanged)
								\@$aref_parameter_pro,	# Start prottest analyses for aa data							-> IN (defined) / OUT (unchanged)
								\@$aref_parameter_mis,	# Replacement code of missing gene sequences 					-> IN (defined) / OUT (unchanged)
								\%$href_outfile_name,	# Outfilename of output option									-> IN (defined) / OUT (unchanged)
					) ;
		} ;
		##############################
		
		
		
		##############################
		## Exclude unchoosen file positions in @$aref_input_files
		my @file_numbers =  split ",", $number_input ;
		my ( %seen_number, @chosen_files ) ;
		
		for( @file_numbers ){ $seen_number{$_}++ } ;
		for ( 0 .. @$aref_input_files-1 ){ push @chosen_files, $aref_input_files->[$_] if $seen_number{$_} }
		
		@$aref_input_files = @chosen_files ;
		
		( %seen_number, @chosen_files, $number_input, $info_text ) = () ;
		##############################
	}
	
	sub input_check{
		
		my $href_hol_seq_of_tax_of_infile		= $_[0] ;	# key: input file name; value ; list of taxa and corresponding sequences				-> IN (undefined) / OUT (defined)
		my $href_taxa_all							= $_[1] ;	# key: taxon name; value : number of occurence in infiles								-> IN (undefined) / OUT (defined)
		my $aref_inputfile_list					= $_[2] ;	# list of input files found in fas, FASTA, phy, or aln format							-> IN (defined) / OUT (unchanged)
		my $href_hoh_info_of_infile_of_type		= $_[3] ;	# key1: infilename, key2: type of info (e.g. 'seqlength'); info type associated value	-> IN (undefined) / OUT (defined)
		my $sref_structure						= $_[4] ;	# defined if structure sequence present in infiles										-> IN (undefined) / OUT (defined)
		
		
		
		##############################
		## If no infile has been found die with an error prompt
		unless ( @$aref_inputfile_list > 0 ){ die "\n\t!FILE-ERROR!: Cannot READ IN infile(s)!\n" }
		##############################
		
		
		
		############################################################ START FILE CHECK
		## READ IN infiles, check for correct format, and make consensus sequence blocks if defined
		for my $file ( sort @$aref_inputfile_list ){
			
			
			##############################
			## counter taxon name, die with an error prompt if taxon name occurs more than once per file
			my %counter_taxa_of_file ;
			##############################
			
			
			
			##############################
			## Print infile READ IN
			print "\n\tREAD IN ", $file;
			##############################
			
			
			
			##############################
			## READ in infiles and assign taxa and corresponding sequences to 
			## %sequence_of_taxon ->  # key: taxon; value sequence of taxon ;
			my %sequence_of_taxon ;
			
			if ( $file =~ /^.*\.aln$/i             ){ &aln2fas ( \$file, \%sequence_of_taxon ) }	# READ IN CLUSTAL files (.aln)
			if ( $file =~ /^.*\.phy$/i             ){ &phy2fas ( \$file, \%sequence_of_taxon ) }	# READ IN PHYLIP  files (.phy)
			if ( $file =~ /^.*\.fas$|^.*\.FASTA$/i ){ &fasEDIT ( \$file, \%sequence_of_taxon ) }	# READ IN FASTA   files (.fas || .FASTA)
			##############################
			
			
			
			##############################
			## Define local variabels for infile chekc
			my	$counter_taxa_of_file ;											# counts specified taxon names in infile. If a taxon name occurs more than once, FASconCAT dies with an error prompt
				$href_hoh_info_of_infile_of_type->{$file}{seqtype} = 'nu' ;	# sequence type predefined as nucleotide data set. If amino acid data unique character states are found, sequence type will be changed to 'aa'
			##############################
			
			
			
			############################################################ START CHECK OF TAXON NAMES AND SEQUENCES
			## Check each taxon and corresponding sequence for...
			## ...equal sequence length like reference sequence length (first sequence in infile)
			## ...correct signs in taxon names (only alphanumeric signs and underscores are allowed in taxon names)
			## ...taxon names with unassigned sequences
			## ...multiple occurence of identic taxon names
			## ...taxon names with associated secondary structure info
			## ...multiple taxa with associated secondary structure info (only one taxon name with secondary structure allowed during concatenation process)
			## ...sequence type (if only one sequence in infile contains amino acid unique caharcter states, infile will be defined as amino acid data set)
			## ...unallowed signs in sequences (only nucleotide and amino acid states as well as missing site, indel, and iupac ambiguity codes are allowed)
			for my $taxon ( sort keys %sequence_of_taxon ){
				
				
				##############################
				## count taxon name, die with an error prompt if taxon name occurs more than once per file
				$counter_taxa_of_file{$taxon}++ ;
				##############################
				
				
				
				##############################
				## if allowed sequence_length of infile is undefined, identify sequence length of given taxon (should be the first infile taxon)
				## and take sequence length as reference sequence length for following sequences of infile
				unless ( $href_hoh_info_of_infile_of_type->{$file}{seqlength} ){ $href_hoh_info_of_infile_of_type->{$file}{seqlength} = length $sequence_of_taxon{$taxon} }
				##############################
				
				
				
				##############################
				## check for...
				## (1) ...unallowed signs in taxon names (only alphanumeric signs and underscores are allowed)
				## (2) ...equal sequence lengths
				## (3) ...missing sequence of identified taxa
				## (4) ...multiple, identical taxon names (not allowed!)
				#die    "\n\t!FILE-ERROR!: Sequence name ", $taxon, " contains unallowed signs in file ", $file, "!\n"	if			$taxon							!~ /^\w+$/ ;				# (1)
				die    "\n\t!FILE-ERROR!: Unequal sequence lengths in file ", $file, "!\n"								if length	$sequence_of_taxon{$taxon}	!= $href_hoh_info_of_infile_of_type->{$file}{seqlength} ;	# (2)
				die    "\n\t!FILE-ERROR!: Sequence missing for taxon ", $taxon, " in file ", $file, "!\n"				unless		$sequence_of_taxon{$taxon}	;						# (3)
				die    "\n\t!FILE-ERROR!: Taxon ", $taxon, " appears multiple times in file ", $file, "!\n"			if			$counter_taxa_of_file{$taxon} > 1 ;					# (4)
				##############################
				
				
				
				############################################################ # START STRUCTURE SEQUENCE IDENTIFICATION & HANDLING
				## Identification of defined secondary structure sequence string...
				## -----------------
				## e.g.: 
				## '(..).((..)..)' or
				## '(--)-((--)--)' or
				## '(..)-((..)..)'
				## -----------------
				## ...and further processing of structure string
				## check for...
				## (1) ...unallowed signs in structure string ( only '(', ')', '.', and '-' allowed )
				## (2) ...multiple taxa with structure strungs (only one taxon with structure string allowed)
				if ( $sequence_of_taxon{$taxon} =~ /.*\(.*\).*/ ){
					
					
					##############################
					## unless taxon name of structure string is defined define given taxon name as structure sequence associated taxon name
					unless ( $href_hoh_info_of_infile_of_type->{$file}{taxstruct} ){ $href_hoh_info_of_infile_of_type->{$file}{taxstruct} = $taxon ; $$sref_structure = $taxon }
					##############################
					
					
					
					##############################
					## substitute structure sign '-' to '.' ('.' -> loop region code)
					$sequence_of_taxon{$taxon}	=~ s/-/./g  ;
					##############################
					
					
					
					##############################
					## Check for correct structure signs ( only '(', ')', '.' are allowed)
					my @strc_elements				=  split "" , $sequence_of_taxon{$taxon} ;
					
					for my $str_sign ( @strc_elements ){
						
						unless ( $str_sign =~ /\(|\)|\./ ){ die "\n\t!FILE-ERROR!: Unallowed signs in structure sequence ", $taxon, " of $file ", $file, "!\n" }	# (1)
					}
					
					@strc_elements = () ;
					##############################
					
					
					
					##############################
					## If given taxon name of associated structure sequence unequal to previous identified structure name, die with an error prompt
					## else print OUT info about identified structure string
					if ( $href_hoh_info_of_infile_of_type->{$file}{taxstruct} ne $taxon ){ die "\n\t!FILE-ERROR!: Multiple structure sequences found: ", $href_hoh_info_of_infile_of_type->{$file}{taxstruct}, ", ", $taxon, "!\n" }	# (2)
					else { print  "\n\t!FILE-INFO!: Structure sequence ", $taxon, " found in file ", $file, "!\n" }
					##############################
				}
				############################################################ # END STRUCTURE SEQUENCE IDENTIFICATION & HANDLING
				
				
				
				############################################################ # START GENE SEQUENCE IDENTIFICATION & HANDLING
				## Identification of gene sequences and further processing
				## substitute lower case character states to upper case states
				## check for...
				## (1) ...unallowed signs in sequences (only nucleotide and amino acid states as well as missing site, indel, and iupac ambiguity codes are allowed)
				## (2) ...sequence type (nucleotide -> 'nu' is predefined. If amino acid unique sequence are found set sequence type of infile to amino acid -> 'aa'
				else{
					
					##############################
					## substitute lower case sequence states to upper case states
					$sequence_of_taxon{$taxon}	=~ s/(\w+)/\U$1/ig ;
					##############################
					
					
					
					##############################
					## Check for unallowed signs in sequences
					my @seq_elements				=  split "" , $sequence_of_taxon{$taxon} ;
					
					for my $seq_sign ( @seq_elements ){ 
						
						unless ( $seq_sign =~ /A|C|G|T|U|-|N|Y|X|R|W|S|K|M|D|V|H|B|Q|E|I|L|F|P|\?/ ){ die "\n\t!FILE-ERROR!: Unknown character found in sequence ", $taxon, " of file ", $file, "!\n" }	# (1)
					}
					
					@seq_elements = () ;
					##############################
					
					
					
					##############################
					## Check for amino acid unique character states
					if ( $sequence_of_taxon{$taxon}  =~ /I|E|L|Q|F|P/i ) { $href_hoh_info_of_infile_of_type->{$file}{seqtype} = 'aa' } # (2)
					##############################
				}
				############################################################ # END GENE SEQUENCE IDENTIFICATION & HANDLING
				
				
				
				##############################
				## Unless consensus sequence has been predefined store...
				## (1) ...taxon and associated sequence (gene/structure string) in hashlist $href_hol_seq_of_tax_of_infile
				## -----------------
				## e.g.: href_hol_seq_of_tax_of_infile->(infile1) = (taxon1, sequence1, taxon2, sequence2, taxon3 ... sequenceN)
				## -----------------
				## (2) ...taxon names of given file as hashkey in %$href_taxa_all; value: number of total occurence in all infiles
				## (3) ...count number of infile taxa
				push @{$href_hol_seq_of_tax_of_infile->{$file}}, ( $taxon, $sequence_of_taxon{$taxon} )	; # (1)
				$href_taxa_all->{$taxon}++																; # (2)
				$href_hoh_info_of_infile_of_type->{$file}{Ntaxa}++										; # (3)
				##############################
			}
			############################################################ END CHECK OF TAXON NAMES AND SEQUENCES
			
			
			
			##############################
			## Print Identified Sequence Type Info
			if 		( $href_hoh_info_of_infile_of_type->{$file}{seqtype} eq 'aa' ){ print "\n\t!FILE-INFO!: Infile ", $file, " identified as AA data set !\n" }
			elsif	( $href_hoh_info_of_infile_of_type->{$file}{seqtype} eq 'nu' ){ print "\n\t!FILE-INFO!: Infile ", $file, " identified as NUC data set !\n" }
			else 	{ die "\n\t!BUG-ERROR!: Cannot assign sequence type of infile ", $file, " in subroutine &input_check!\n\tPlease, report BUG to system developer!" }
			##############################
		}
		############################################################ END FILE CHECK
		
		
		############################################################ START SUB-SUBROUTINES WITHIN SUBROUTINE &input_check
		# READ IN CLUSTAL FORMATED FILES
		sub aln2fas{
			
			my $sref_file_aln		= $_[0] ;	# infile name							-> IN (defined) / OUT (unchanged)
			my $href_seq_of_tax	= $_[1] ;	# key: taxon; value sequence of taxon	-> IN (undefined) / OUT (defined)
			
			
			
			####################################
			## Untie line feeds
			&tie_linefeeds( \$$sref_file_aln ) ;
			####################################
			
			
			
			my (
				@tax_seq_split,
				@taxa,
				$clustal,
				%sequenzen_alle,
				@inputfile
			);
			
			
			####################################
			## READ IN clustal file
			open INaln, $$sref_file_aln || die "\n\tFILE-ERROR: Cannot READ IN ", $$sref_file_aln, "!\n" ;
			
			while (my $line = <INaln>){
				
				chomp  $line ;
				if   ( $line =~ /^CLUSTAL/i ){ $clustal=1 }
				push   @inputfile, "$line\n" ;
			}
			close INaln;
			####################################
			
			if ($clustal == 1) { splice (@inputfile, 0, 2) } else { die "\n\tFILE-ERROR: Unknown clustal format found in file ", $$sref_file_aln, "!\n" }
			
			for (@inputfile){ if ( /^\W/ ) { s/.*/:/g ; s/\n// } }
			
			my	$string_inputfile	=  join ""		, @inputfile ;
				$string_inputfile	=~ s/:+/:/g		, $string_inputfile ;
			my @seq_parts			=  split ":"	, $string_inputfile ;
			
			for ( @seq_parts ){ 
				
				my @tax_seq = split "\n", $_ ;
				
				for ( @tax_seq ){ s/ +/:/ ; push @tax_seq_split, (split ":", $_) }
				
				my %sequenzen  = @tax_seq_split ;
				@tax_seq_split = () ;
				@taxa          = keys %sequenzen ;
				
				for my $taxon( @taxa ){ push(@{$sequenzen_alle{$taxon}}, $sequenzen{$taxon}) }
			}
			
			for my $taxon_aln(@taxa){
				
				my   $aref_a = exists($sequenzen_alle{$taxon_aln}) ? \@{$sequenzen_alle{$taxon_aln}} :() ;
				my   $sequenz_final = join "", @$aref_a ;
				$href_seq_of_tax->{$taxon_aln} = $sequenz_final ;
			}
		}
		
		# READ IN PHYLIP FORMATED FILES
		sub phy2fas{
			
			my $sref_phy_file		= $_[0] ;	# infile name							-> IN (defined) / OUT (unchanged)
			my $href_seq_of_tax	= $_[1] ;	# key: taxon; value sequence of taxon	-> IN (undefined) / OUT (defined)
			
			
			
			####################################
			## Untie line feeds
			&tie_linefeeds( \$$sref_phy_file ) ;
			####################################
			
			
			
			########################################################################
			## READ IN phylip interleaved and non-interleaved formated files
			## store single taxa and their associated sequence in has $href_seq_of_tax
			open INphy , $$sref_phy_file or die "\n\t!FILE-ERROR!: Cannot READ IN file ", $$sref_phy_file, " !\n" ;
			chomp (my @all_lines_phy = <INphy>) and close INphy ; 
			
			
			
			####################################
			## extract the first line (infoline) and determine the number of taxa ($info_line[0])
			( my $infoline	= shift @all_lines_phy ) =~ s/\s+/ / ;
			my	@infos_line	=  split " ", $infoline ;
			####################################
			
			
			
			####################################
			## phylip files can be in interleaved and non-interleaved format
			## interleaved:	tax1 ACGT...	# first part of all taxa
			## 				tax2 ACGT...
			##								#''space line''
			##				tax1 CCCC...	# Second part of all taxa
			##				tax2 GGGG...
			## to concatenate sequence parts correctly store single lines in separate hashkeys
			## until number of taxa is reached ($c equal $infos_line[0]). afterwards remove the following spaceline and concatenate
			## next lines to their corresponding taxon sequences inferred from the first round and so on...
			## If phylip file is in non-interleaved format, the while lopp stops automatically after the first foreach loop
			my 	%seq_phy = () ;
			while ( @all_lines_phy ){
				
				for ( my $c=1; $c<=$infos_line[0]; $c++ ){ my $seq_line_phy = shift @all_lines_phy ; push ( @{$seq_phy{$c}} , $seq_line_phy ) }
				shift @all_lines_phy ;
			}
			####################################
			
			
			
			####################################
			## join single sequence parts of each taxon (if interleaved there are multiple key values)
			## taxonnames are in the same line as sequenceinformation (separated by one or multiple whitespaces), therefore
			## substitute multiple whitespaces into one whitespace, join all sequence parts to one string (only important if interleaved format)
			## split taxonnames from sequenceinformation and store both in the hashreference %href_seq_of_tax (key: taxon; value: complete sequence)
			my %seen_phylip_taxon ;
			for my $line_c ( sort {$a<=>$b} keys %seq_phy ){ 
				
				my	@seq_single_parts				=	exists($seq_phy{$line_c}) ? @{$seq_phy{$line_c}} :( ) ;
					$seq_single_parts[0]			=~	s/\s+/ / ;
				my	$seq_complete					=	join	""		, @seq_single_parts ;
					@seq_single_parts				=	split	" "		, $seq_complete ;
					
					if ( @seq_single_parts > 2 ){ die "\n\tFILE-ERROR!: At least one taxon name contains blanks in phylip file ", $$sref_phy_file, "! Not allowed for phylip files!\n\n"} 
				my	$taxon							=	shift @seq_single_parts ;
				
				unless ( $seen_phylip_taxon{$taxon} ){ $seen_phylip_taxon{$taxon}++ } else{ die "\n\t!FILE-ERROR!: Taxon ", $taxon, " appears multiple times in file ", $$sref_phy_file, " !\n" }
				
					$seq_complete					=	join	""		, @seq_single_parts ;
					$href_seq_of_tax->{$taxon}	=	$seq_complete ;
			}
			########################################################################
		}
		
		# READ IN FASTA FORMATED FILES
		sub fasEDIT{
			
			my $sref_file_fas		= $_[0] ; # infile name							-> IN (defined) / OUT (unchanged)
			my $href_seq_of_tax	= $_[1] ; # key: taxon; value sequence of taxon	-> IN (undefined) / OUT (defined)
			
			
			
			####################################
			## Untie line feeds
			&tie_linefeeds( \$$sref_file_fas ) ;
			####################################
			
			
			
			####################################
			## READ IN fasta interleaved and non-interleaved formated files
			## store single taxa and their associated sequence in has $href_sequence_of_taxon
			my 	(
					$taxon,			# taxon name
					%seen_fasta_taxon,	# key: taxon name; value counter of taxon occurence in infile
			) ;
			
			open INfas, $$sref_file_fas || die "\n\tFILE-ERROR: Cannot READ IN ", $$sref_file_fas, "!\n" ;
			
			while ( my $line = <INfas> ){
				
				chomp	$line ;
				
				if ( $line =~ /^\>/ ){
					
					( $taxon = $line ) =~ s/^\>// ;
					
					unless ($seen_fasta_taxon{$taxon}){ $seen_fasta_taxon{$taxon}++ }
					else{ die "\n\t!FILE-ERROR!: Taxon ", $taxon, " appears multiple times in fasta file ", $$sref_file_fas, " !\n" }
				}
				 
				else { 
					
					if ( $seen_fasta_taxon{$taxon} ){ $href_seq_of_tax->{$taxon} .= $line }
					else { die "\n\t!FILE-ERROR!: Cannot assign first sequence in fasta file ", $$sref_file_fas, " !\n" }
				}
			}
			
			close INfas
			####################################
		}
		
		# CHANGE LINEFEEDS OF INPUT LINES TO LINUX STANDARD (LF)
		sub tie_linefeeds{
			
			my $sref_tie_file = $_[0] ;
			
			# Untei linefeeds
			TIE:
			(tie ( my @data, 'Tie::File', $$sref_tie_file )) ;
			
			print  "\n\t!FILE-ERROR!: ", $$sref_tie_file, " is empty!\n" and next READING if 0 == @data ;
			
			map { s/\r\n/\n/g } @data ;
			map { s/\r/\n/g   } @data ;
			
			untie @data ;
		}
		############################################################ END SUB-SUBROUTINES WITHIN SUBROUTINE &input_check
	} 
	
	sub seq_translation{
		
		my $href_hol_seq_of_tax_of_infile		= $_[0] ;	# key: input file name; value ; list of taxa and corresponding sequences				-> IN (defined) / OUT (changed)
		my $href_hoh_info_of_infile_of_type		= $_[1] ;	# key1: infilename, key2: type of info (e.g. 'seqlength'); info type associated value	-> IN (defined) / OUT (changed)
		my $sref_tra_parameter					= $_[2] ;	# defined ry coding parameter option													-> IN (defined) / OUT (unchanged)
		
		
		##############################
		# Definition of forward and reward triplet code
		# forward:	-> key: amino acid state; value: nucleotide triplet
		# reward	-> key: nucleotide triplet; value: amino acid code
		my %triplet_code = (
							'forward'	=>
							{
									'F' => 'TTY',
									'A' => 'GCN',
									'R' => 'MGR',
									'N' => 'AAY',
									'D' => 'GAY',
									'C' => 'TGY',
									'Q' => 'CAR',
									'E' => 'GAR',
									'G'	=> 'GGN',
									'H'	=> 'CAY',
									'I'	=> 'ATH',
									'M'	=> 'ATG',
									'K'	=> 'AAR',
									'P'	=> 'CCN',
									'L'	=> 'YTR',
									'S'	=> 'WSY',
									'T'	=> 'ACN',
									'W'	=> 'TGG',
									'Y'	=> 'TAY',
									'V'	=> 'GTN'
							},
							
							'reward'	=>
							{
									'TTT' => 'F',	# Phenylalanin
									'TTC' => 'F',
									'TTY' => 'F',
									
									'TTA' => 'L',	# Leucine
									'TTG' => 'L',
									'TTR' => 'L',
									'CTT' => 'L',
									'CTC' => 'L',
									'CTA' => 'L',
									'CTG' => 'L',
									'CTM' => 'L',
									'CTR' => 'L',
									'CTW' => 'L',
									'CTS' => 'L',
									'CTY' => 'L',
									'CTK' => 'L',
									'CTV' => 'L',
									'CTH' => 'L',
									'CTD' => 'L',
									'CTB' => 'L',
									'CTN' => 'L',
									'YTT' => 'L',
									'YTC' => 'L',
									'YTA' => 'L',
									'YTG' => 'L',
									'YTR' => 'L',
									
									'ATT' => 'I',	# Isoleucine
									'ATC' => 'I',
									'ATA' => 'I',
									'ATH' => 'I',
									
									'ATG' => 'M',	# Methionine
									
									'GTT' => 'V',	# Valine
									'GTC' => 'V',
									'GTA' => 'V',
									'GTG' => 'V',
									'GTM' => 'V',
									'GTR' => 'V',
									'GTW' => 'V',
									'GTS' => 'V',
									'GTY' => 'V',
									'GTK' => 'V',
									'GTV' => 'V',
									'GTH' => 'V',
									'GTD' => 'V',
									'GTB' => 'V',
									'GTN' => 'V',
									
									'TCT' => 'S',	# Serine
									'TCC' => 'S',
									'TCA' => 'S',
									'TCG' => 'S',
									'TCM' => 'S',
									'TCR' => 'S',
									'TCW' => 'S',
									'TCS' => 'S',
									'TCY' => 'S',
									'TCK' => 'S',
									'TCV' => 'S',
									'TCH' => 'S',
									'TCD' => 'S',
									'TCB' => 'S',
									'TCN' => 'S',
									'AGT' => 'S',
									'AGC' => 'S',
									'AGY' => 'S',
									'WSY' => 'S',
									
									'CCT' => 'P',	# Proline
									'CCC' => 'P',
									'CCA' => 'P',
									'CCG' => 'P',
									'CCM' => 'P',
									'CCR' => 'P',
									'CCW' => 'P',
									'CCS' => 'P',
									'CCY' => 'P',
									'CCK' => 'P',
									'CCV' => 'P',
									'CCH' => 'P',
									'CCD' => 'P',
									'CCB' => 'P',
									'CCN' => 'P',
									
									'ACT' => 'T',	# Threonine
									'ACC' => 'T',
									'ACA' => 'T',
									'ACG' => 'T',
									'ACM' => 'T',
									'ACR' => 'T',
									'ACW' => 'T',
									'ACS' => 'T',
									'ACY' => 'T',
									'ACK' => 'T',
									'ACV' => 'T',
									'ACH' => 'T',
									'ACD' => 'T',
									'ACB' => 'T',
									'ACN' => 'T',
									
									'GCT' => 'A',	# Alanine
									'GCC' => 'A',
									'GCA' => 'A',
									'GCG' => 'A',
									'GCM' => 'A',
									'GCR' => 'A',
									'GCW' => 'A',
									'GCS' => 'A',
									'GCY' => 'A',
									'GCK' => 'A',
									'GCV' => 'A',
									'GCH' => 'A',
									'GCD' => 'A',
									'GCB' => 'A',
									'GCN' => 'A',
									
									'TAT' => 'Y',	# Tyrosine
									'TAC' => 'Y',
									'TAY' => 'Y',
									
									'CAT' => 'H',	# Histidine
									'CAC' => 'H',
									'CAY' => 'H',
									
									'CAA' => 'Q',	# Glutamine
									'CAG' => 'Q',
									'CAR' => 'Q',
									
									'AAT' => 'N',	# Asparagine
									'AAC' => 'N',
									'AAY' => 'N',
									
									'AAA' => 'K',	# Lysine
									'AAG' => 'K',
									'AAR' => 'K',
									
									'GAT' => 'D',	# Aspartic acid
									'GAC' => 'D',
									'GAY' => 'D',
									
									'GAA' => 'E',	# Glutamic acid
									'GAG' => 'E',
									'GAR' => 'E',
									
									'TGT' => 'C',	# Cysteine
									'TGC' => 'C',
									'TGY' => 'C',
									
									'TGG' => 'W',	# Tryptophan
									
									'CGT' => 'R',	# Arginine
									'CGC' => 'R',
									'CGA' => 'R',
									'CGG' => 'R',
									'CGM' => 'R',
									'CGR' => 'R',
									'CGW' => 'R',
									'CGS' => 'R',
									'CGY' => 'R',
									'CGK' => 'R',
									'CGV' => 'R',
									'CGH' => 'R',
									'CGD' => 'R',
									'CGB' => 'R',
									'CGN' => 'R',
									'AGA' => 'R',
									'AGG' => 'R',
									'AGR' => 'R',
									'MGA' => 'R',
									'MGG' => 'R',
									'MGR' => 'R',
									
									'GGT' => 'G',	# Glycine
									'GGC' => 'G',
									'GGA' => 'G',
									'GGG' => 'G',
									'GGM' => 'G',
									'GGR' => 'G',
									'GGW' => 'G',
									'GGS' => 'G',
									'GGY' => 'G',
									'GGK' => 'G',
									'GGV' => 'G',
									'GGH' => 'G',
									'GGD' => 'G',
									'GGB' => 'G',
									'GGN' => 'G',
							},
		) ;
		##############################
		
		
		
		############################################################ START SEQUENCE TRANSLATION OF INFILES
		for my $infile ( sort {$a<=>$b} keys %$href_hol_seq_of_tax_of_infile ){
			
			
			##############################
			## Store identified sequences and associated taxon names of given infile in...
			## ...%sequence_of_taxon -> key: original taxon name; value: associated sequence
			my $aref_list_taxon_sequence_of_file	= exists( $href_hol_seq_of_tax_of_infile->{$infile}) ? \@{$href_hol_seq_of_tax_of_infile->{$infile}} :( ) ;
			my %sequence_of_taxon						= @$aref_list_taxon_sequence_of_file ;
			##############################
			
			
			
			############################################################ START SEQUENCE TRANSLATION OF NUCLEOTIDE INFILE TO AMINO ACID DATA (IF DEFINED)
			if ( ( $href_hoh_info_of_infile_of_type->{$infile}{seqtype} eq 'nu' ) && ( $$sref_tra_parameter eq 'NUC to AA' ) ){
				
				
				##############################
				## Infiles with secondary structure info are not translated to AA data
				unless ( $href_hoh_info_of_infile_of_type->{$infile}{taxstruct} ){
					
					
					##############################
					## Print file handling info
					print "\n\tTranslating NUC sequences of ", $infile, " to AA sequences..." ;
					##############################
					
					
					
					##############################
					## Print warning if sequence lengths are not a multiple of three
					if ( ( $href_hoh_info_of_infile_of_type->{$infile}{seqlength} % 3 ) != 0 ){ print "\n\tINFO: Sequence lengths of infile ", $infile, " not a multiple of 3!\n" }
					##############################
					
					
					
					##############################
					## For each taxon start triplett substitution to amino acid state
					## As long as sequence states are given in nucleotide sequence...
					## (1) ...extract substrings of 3 states and
					## (2) ...assign, if defined triplet, corresponding amino acid state to new sequence; else assign '?'
					## (3) ...store translated sequence as hashvalue of given taxon
					## (4) ...store new sequence length as hashvalue of give infile
					my $new_seq_length ;
					for my $taxon ( keys %sequence_of_taxon ){
						
						my $translated_seq ;
						for ( my $k=0; $k<=$href_hoh_info_of_infile_of_type->{$infile}{seqlength}-1; $k+=3 ){
							
							# (1)
							my $triplet_nuc = substr($sequence_of_taxon{$taxon}, $k, 3) ;
							##
							
							# (2)
							if ( $triplet_code{reward}{$triplet_nuc} ){ $translated_seq .= $triplet_code{reward}{$triplet_nuc} } else{ $translated_seq .= '?' }
							##
							
							# test print
							#print " ", $triplet_nuc, "\t", $translated_seq, "\n"
							##
						}
						
						# (3)
						$sequence_of_taxon{$taxon} = $translated_seq ;
						##
						
						# (4)
						unless ( $new_seq_length ){ $new_seq_length = length $translated_seq }
						##
					}
					##############################
					
					
					
					##############################
					## Assign new sequence type and sequence length of infile
					$href_hoh_info_of_infile_of_type->{$infile}{seqtype}	= 'aa' ;
					$href_hoh_info_of_infile_of_type->{$infile}{seqlength}	= $new_seq_length ;
					##############################
				}
				else{ print "\n\tCannot translate NUC sequences of ", $infile, " to AA sequences due to identified structure string info!" ; }
			}
			############################################################ END SEQUENCE TRANSLATION OF NUCLEOTIDE INFILE TO AMINO ACID DATA (IF DEFINED)
			
			
			
			############################################################ START SEQUENCE TRANSLATION OF AMINO ACID DATA INFILE TO NUCLEOTIDE (IF DEFINED)
			elsif ( ( $href_hoh_info_of_infile_of_type->{$infile}{seqtype} eq 'aa' ) && ( $$sref_tra_parameter eq 'AA to NUC' ) ){
				
				
				##############################
				## Print file handling info
				print "\n\tTranslating AA sequences of ", $infile , " to NUC sequences...";
				##############################
				
				
				
				##############################
				## for each taxon start amino acid substitution to to nucleotide triplet
				## As long as sequence states are given in nucleotide sequence...
				## (1) ...extract substrings of 1 state and
				## (2) ...substitute state, if defined amino acid state, to corresponding nucleotide triplet; else assign '???' to new sequence
				## (3) ...store translated sequence as hashvalue of given taxon
				for my $taxon ( keys %sequence_of_taxon ){
					
					my $translated_seq ;
					for ( my $k=0; $k<=$href_hoh_info_of_infile_of_type->{$infile}{seqlength}-1; $k+=1 ){
						
						# (1)
						my $aa_site = substr($sequence_of_taxon{$taxon}, $k, 1) ;
						##
						
						# (2)
						if ( $triplet_code{forward}{$aa_site} ){ $translated_seq .= $triplet_code{forward}{$aa_site} } else{ $translated_seq .= '???' }
						##
						
						# test print
						#print " ", $aa_site, "\t", $translated_seq, "\n"
						##
					}
					
					# (3)
					$sequence_of_taxon{$taxon} = $translated_seq
					##
				}
				##############################
				
				
				
				##############################
				## Assign new sequence type of infile
				$href_hoh_info_of_infile_of_type->{$infile}{seqtype} = 'nu'
				##############################
			}
			############################################################ START SEQUENCE TRANSLATION OF AMINO ACID DATA INFILE TO NUCLEOTIDE (IF DEFINED)
			
			
			
			##############################
			## assign recoded sequences to %$href_hol_seq_of_tax_of_infile
			$href_hol_seq_of_tax_of_infile->{$infile} = () ;
			
			for my $taxon ( sort keys %sequence_of_taxon ){ push @{$href_hol_seq_of_tax_of_infile->{$infile}}, ( $taxon, $sequence_of_taxon{$taxon} ) }
			##############################
		}
		############################################################ END SEQUENCE TRANSLATION OF INFILES
	}
	
	sub ry_coding{
		
		my $href_hol_seq_of_tax_of_infile		= $_[0] ;	# key: input file name; value ; list of taxa and corresponding sequences				-> IN (defined) / OUT (changed)
		my $href_hoh_info_of_infile_of_type		= $_[1] ;	# key1: infilename, key2: type of info (e.g. 'seqlength'); info type associated value	-> IN (defined) / OUT (unchanged)
		my $sref_ry_parameter					= $_[2] ;	# defined ry coding parameter option													-> IN (defined) / OUT (unchanged)
		
		
		
		my %ry_code_of_state = (
									'A' => 'R',
									'G' => 'R',
									'C' => 'Y',
									'T' => 'Y',
									'U' => 'Y',
									'?' => '?',
									'-' => '-',
									'M'	=> 'M',
									'R'	=> 'R',
									'W'	=> 'W',
									'S'	=> 'S',
									'Y'	=> 'Y',
									'K'	=> 'K',
									'V'	=> 'V',
									'H'	=> 'H',
									'D'	=> 'D',
									'B'	=> 'B',
									'N'	=> 'N'
		) ;
		
		
		############################################################ START RY CODING OF EACH INFILE
		for my $infile ( sort {$a<=>$b} keys %$href_hol_seq_of_tax_of_infile ){
			
			
			
			############################################################ START HANDLING OF NUC DATA
			if ( $href_hoh_info_of_infile_of_type->{$infile}{seqtype} eq 'nu' ){
				
				
				##############################
				## Print file handling info
				print "\n\tRY coding of NUC sequences in ", $infile, "..." ;
				##############################
				
				
				
				##############################
				## Store identified sequences and associated taxon names of given infile in...
				## ...%sequence_of_taxon -> key: original taxon name; value: associated sequence
				my $aref_list_taxon_sequence_of_file	= exists( $href_hol_seq_of_tax_of_infile->{$infile}) ? \@{$href_hol_seq_of_tax_of_infile->{$infile}} :( ) ;
				my %sequence_of_taxon						= @$aref_list_taxon_sequence_of_file ;
				my %RY_code_of_taxon;
				##############################
				
				
				
				##############################
				## Recode taxon sequences
				for my $taxon ( sort keys %sequence_of_taxon ){
					
					
					############################################################ START RY CODING OF EACH NUC TAXON SEQ OF INFILE
					unless ( $href_hoh_info_of_infile_of_type->{$infile}{taxstruct} eq $taxon ){
						
						my @states = split "", $sequence_of_taxon{$taxon} ;
						my $recoded_sequence ;
						
						
						############################################################ START RECODING 3rd SEQ POSITIONS
						if ( $$sref_ry_parameter eq '3rd' ){
							
							while ( @states ){
								
								my $first	= shift @states ;
								my $second	= shift @states ;
								my $third	= shift @states ;
								
								$recoded_sequence .= $first.$second.$ry_code_of_state{$third}
							}
						}
						############################################################ END RECODING 3rd SEQ POSITIONS
						
						
						
						############################################################ START RECODING ALL SEQ POSITIONS
						elsif ( $$sref_ry_parameter eq 'All' ){
							
							for my $state ( @states ){ $recoded_sequence .= $ry_code_of_state{$state} }
							@states = () ;
						}
						############################################################ END RECODING ALL SEQ POSITIONS
						
						
						
						##############################
						## if defined ry coding option cannot be assigned, print bug report
						else{ die "\n\tBUG-ERROR: Cannot assign RY coding option ", $$sref_ry_parameter, "!\n\tPlease, report BUG to system developer!" } 
						##############################
						
						
						
						##############################
						## assign new sequence of taxon
						$RY_code_of_taxon{$taxon} = $recoded_sequence ;
						##############################
					}
					############################################################ END RY CODING OF EACH NUC TAXON SEQ OF INFILE
				}
				##############################
				
				
				
				##############################
				## assign recoded sequences to %$href_hol_seq_of_tax_of_infile
				$href_hol_seq_of_tax_of_infile->{$infile}  = () ;
				
				for my $taxon ( sort keys %RY_code_of_taxon ){ push @{$href_hol_seq_of_tax_of_infile->{$infile}}, ( $taxon, $RY_code_of_taxon{$taxon} ) }
				
				%RY_code_of_taxon = () ;
				%sequence_of_taxon = () ;
				##############################
			}
			############################################################ END HANDLING OF NUC DATA
			
			
			
			##############################
			## KEEP AA DATA unchanged
			else{ print "\n\t", $infile, " consists AA data, sequences not RY coded!" }
			##############################
		}
		############################################################ END RY CODING OF EACH INFILE
	}
	
	sub reject_third_nuc_position{
		
		my $href_hol_seq_of_tax_of_infile		= $_[0] ;	# key: input file name; value ; list of taxa and corresponding sequences				-> IN (defined) / OUT (changed)
		my $href_hoh_info_of_infile_of_type		= $_[1] ;	# key1: infilename, key2: type of info (e.g. 'seqlength'); info type associated value	-> IN (defined) / OUT (changed)
		
		
		############################################################ START EXCLUSION OF THIRD SEQ PSITION IN EACH INFILE
		for my $infile ( sort {$a<=>$b} keys %$href_hol_seq_of_tax_of_infile ){
			
			
			##############################
			## Store identified sequences and associated taxon names of given infile in...
			## ...%sequence_of_taxon -> key: original taxon name; value: associated sequence
			my $aref_list_taxon_sequence_of_file	= exists( $href_hol_seq_of_tax_of_infile->{$infile}) ? \@{$href_hol_seq_of_tax_of_infile->{$infile}} :( ) ;
			my %sequence_of_taxon						= @$aref_list_taxon_sequence_of_file ;
			##############################
			
			
			
			############################################################ START EXCLUSION OF THIRD NUC SEQUENCE POSITION
			## If sequence type = nucleotide data, remove third sequence position
			## Print Info if sequence lengths of infile are not a multiple of three
			## Store new sequence length in $href_hoh_info_of_infile_of_type->{$infile}{seqlength}
			## Store new sequences in %$href_hol_seq_of_tax_of_infile
			if ( $href_hoh_info_of_infile_of_type->{$infile}{seqtype} eq 'nu' ){
				
				
				##############################
				## Print warning if sequence lengths are not a multiple of three
				if ( ( $href_hoh_info_of_infile_of_type->{$infile}{seqlength} % 3 ) != 0 ){ print "\n\tINFO: Sequence lengths of infile ", $infile, " not a multiple of 3!" }
				##############################
				
				
				
				##############################
				## Print file handling info
				print "\n\tExcluding each 3rd NUC site of ", $infile, "..." ;
				##############################
				
				
				
				##############################
				## Reduce sequences of each taxon
				for my $taxon ( sort keys %sequence_of_taxon ){
					
					
					##############################
					## store reduced sequence in...
					my $reduced_sequence ;
					##############################
					
					
					
					############################## NUC DATA
					## Reduction of nucleotide data
					unless ( $href_hoh_info_of_infile_of_type->{$infile}{taxstruct} eq $taxon ){
						
						for ( my $j=0; $j<= $href_hoh_info_of_infile_of_type->{$infile}{seqlength}; $j+=3 ){ $reduced_sequence .= substr($sequence_of_taxon{$taxon}, $j, 2) }
						
						## test print
						#print "\n", $taxon, "\t", $reduced_sequence, "\n"; exit;
						####
					}
					##############################
					
					
					
					############################## STRUC DATA
					## Reduction of structure sequence data...
					## (1) ...identify paired stem sequence positions
					## (2) ...substitute stem position to loop '.' if counterpart is a multiple of three
					## (3) ...exclude third sequence positions
					else{
						
						my @states = split "", $sequence_of_taxon{$taxon} ;
						
						## (1)
						my ( @open_brackets, @pairs ) ;
						for ( 0 .. @states-1 ){
							
							if 		( $states[$_] =~ /\(/ ){ push @open_brackets, $_ }
							elsif 	( $states[$_] =~ /\)/ ){ my $open_pos = pop @open_brackets; push @pairs, $open_pos.':'.$_ }
						}
						####
						
						
						
						## (2)
						for ( @pairs ){
							
							my @positions = split ":", $_ ;
							if 		( $positions[0] % 3 ){ $states[$positions[1]] = "." }
							elsif 	( $positions[1] % 3 ){ $states[$positions[0]] = "." }
						}
						####
						
						
						
						## (3)
						my $orig_seq = join "", @states ;
						for ( my $j=0; $j<= $href_hoh_info_of_infile_of_type->{$infile}{seqlength}; $j+=3 ){ $reduced_sequence .= substr($orig_seq, $j, 2) }
						####
					}
					##############################
					
					
					
					##############################
					## assign new sequence of taxon
					$sequence_of_taxon{$taxon} = $reduced_sequence ;
					##############################
				}
				##############################
				
				
				
				##############################
				## assign reduced sequences to %$href_hol_seq_of_tax_of_infile
				## assign new sequence lengths of infile
				$href_hol_seq_of_tax_of_infile->{$infile} = () ;
				
				for my $taxon ( sort keys %sequence_of_taxon ){ push @{$href_hol_seq_of_tax_of_infile->{$infile}}, ( $taxon, $sequence_of_taxon{$taxon} ); $href_hoh_info_of_infile_of_type->{$infile}{seqlength} = length $sequence_of_taxon{$taxon} }
				##############################
			}
			############################################################ END EXCLUSION OF THIRD NUC SEQUENCE POSITION
			
			
			
			##############################
			## KEEP THIRD SEQUENCE POSITION OF AA DATA
			else{ print "\n\tInfile ", $infile, " consists AA data, keep 3rd site position!\n" }
			##############################
		}
		############################################################ END EXCLUSION OF THIRD SEQ PSITION IN EACH INFILE
	}
	
	sub concatenate{
		
		my $href_hol_seq_of_tax_of_infile		= $_[0] ;	# key: input file name; value ; list of taxa and corresponding sequences				-> IN (defined) / OUT (changed)
		my $href_hoh_info_of_infile_of_type		= $_[1] ;	# key1: infilename, key2: type of info (e.g. 'seqlength'); info type associated value	-> IN (undefined) / OUT (defined)
		my $href_name_list						= $_[2] ;	# key: taxon name; value : number of occurence in infiles								-> IN (defined) / OUT (unchanged)
		my $sref_structure						= $_[3] ;	# defined if structure sequence present in infiles										-> IN (defined) / OUT (unchanged)
		my $sref_fill_code						= $_[4] ;	# defined replacement parameter (MISSING or INDEL)										-> IN (defined) / OUT (unchanged)
		
		
		##############################
		## Define supermatrix sequence type as nucleotide data type
		## unless single infile sequence has been identified as amino acid data type
		$href_hoh_info_of_infile_of_type->{supermatrix}{seqtype} = 'nu' ;
		##############################
		
		
		
		########################################################################################## START CONCATENATION OF EACH INFILE & SEQUENCE INFO EXTRACTION
		## For each infile...
		## ...determine start and endposition in supermatrix file
		## ...extract additional sequence info
		## ...concatenate single infile sequences due to corresponding taxon names
		my %supermatrix_of_taxon ;
		$href_hoh_info_of_infile_of_type->{supermatrix}{seqstart} = 1 ;
		$href_hoh_info_of_infile_of_type->{supermatrix}{seqend	} = 1 ;
		
		my $counter_bp = 1 ;
		for my $infile ( sort keys %$href_hol_seq_of_tax_of_infile ){
			
			
			
			##############################
			## KEEP AA DATA unchanged
			print "\n\tConcatenate ", $infile, "..." ;
			##############################
			
			
			
			##############################
			## Store structure sequence name in supermatrix infile
			if ( $href_hoh_info_of_infile_of_type->{$infile}{taxstruct}	){
				
				$href_hoh_info_of_infile_of_type->{supermatrix}{taxstruct} = $href_hoh_info_of_infile_of_type->{$infile}{taxstruct}
			}
			##############################
			
			
			
			##############################
			## if single infile amino acid data -> set supermatrix type to amino acid data
			if ( $href_hoh_info_of_infile_of_type->{$infile}{seqtype} eq 'aa' ){ $href_hoh_info_of_infile_of_type->{supermatrix}{seqtype} = 'aa' }
			##############################
			
			
			
			##############################
			## Store identified sequences and associated taxon names of given infile in...
			## ...%sequence_of_taxon -> key: original taxon name; value: associated sequence
			my $aref_list_taxon_sequence_of_file	= exists( $href_hol_seq_of_tax_of_infile->{$infile}) ? \@{$href_hol_seq_of_tax_of_infile->{$infile}} :( ) ;
			my %sequence_of_taxon					= @$aref_list_taxon_sequence_of_file ;
			##############################
			
			
			
			############################################################ START DETERMINATION OF START & ENDPOSITION OF INFILE SEQUENCES IN SUPERMATRIX
			## Set supermatrix start position to $counter_bp
			## If sequence length > 2 base pairs, set supermatrix endposition = start position + sequence length of gene partition -1 
			## else set supermatrix endposition to $counter_bp
			$href_hoh_info_of_infile_of_type->{$infile}{seqstart}		= $href_hoh_info_of_infile_of_type->{supermatrix}{seqend} ;
			
			if ( $href_hoh_info_of_infile_of_type->{$infile}{seqlength} == 1 ){
				
				$href_hoh_info_of_infile_of_type->{$infile}{seqend}		= $href_hoh_info_of_infile_of_type->{supermatrix}{seqend}
			}
			else{
				
				$href_hoh_info_of_infile_of_type->{$infile}{seqend}		= $href_hoh_info_of_infile_of_type->{$infile}{seqstart} + $href_hoh_info_of_infile_of_type->{$infile}{seqlength} - 1 ;
				$href_hoh_info_of_infile_of_type->{supermatrix}{seqend}	= $href_hoh_info_of_infile_of_type->{$infile}{seqstart} + $href_hoh_info_of_infile_of_type->{$infile}{seqlength}
			}
			############################################################ END DETERMINATION OF START & ENDPOSITION OF INFILE SEQUENCES IN SUPERMATRIX
			
			
			
			############################################################ START CONCATENATION PROCESS
			## Sequence concatenation process
			## if sequence of taxon name found in infile, elongate supermatrix sequence of that taxon
			## elsif sequence of taxon not given in infile, elongate taxon sequence with 
			## ambiguity characters...
			## (1) ... dots if taxon sequence includes structure sequence
			## (2) ... indels ('-') if specified via -x option, otherwise with
			## (3) ... N's if infile sequences consist of nucleotide sequence
			## (4) ... X's if infile sequences consist of amino acid sequence
			##
			## Store taxon with sequences in infile in $href_hol_seq_of_tax_of_infile->{$infile}{Nseqconcat}
			## Store missing taxon names of infile in
			my (
					$string_missing_taxon_dots,		#
					$string_missing_taxon_states,	#
			) ;
			
			if 		( $$sref_structure														)	{ $string_missing_taxon_dots		= ( "." x $href_hoh_info_of_infile_of_type->{$infile}{seqlength} ) }
			if 		( $$sref_fill_code										eq 'Indel' 		)	{ $string_missing_taxon_states		= ( "-" x $href_hoh_info_of_infile_of_type->{$infile}{seqlength} ) }
			elsif	( $href_hoh_info_of_infile_of_type->{$infile}{seqtype}	eq 'nu' 		)	{ $string_missing_taxon_states		= ( "N" x $href_hoh_info_of_infile_of_type->{$infile}{seqlength} ) }
			else 																				{ $string_missing_taxon_states		= ( "X" x $href_hoh_info_of_infile_of_type->{$infile}{seqlength} ) }
			
			for my $taxon ( sort {$a<=>$b} keys %$href_name_list ){
				
				$href_hoh_info_of_infile_of_type->{supermatrix}{taxconcat}{$taxon}++ ;
				
				if 		( $sequence_of_taxon{$taxon} 											){ $supermatrix_of_taxon{$taxon} .= $sequence_of_taxon{$taxon} 		; $href_hoh_info_of_infile_of_type->{supermatrix}{Nseqconcat  }{$taxon}++	; $href_hoh_info_of_infile_of_type->{$infile}{Ntaxa_in}++		}
				elsif	( $taxon	eq $$sref_structure											){ $supermatrix_of_taxon{$taxon} .= $string_missing_taxon_dots		; $href_hoh_info_of_infile_of_type->{supermatrix}{Nsequnconcat}{$taxon}++	; $href_hoh_info_of_infile_of_type->{$infile}{Ntaxa_mis}++	} #push @{$href_hol_seq_of_tax_of_infile->{$infile}}, ( $taxon, $string_missing_taxon_dots	) }	# (1)
				elsif	( 'nu'		eq $href_hoh_info_of_infile_of_type->{$infile}{seqtype  }	){ $supermatrix_of_taxon{$taxon} .= $string_missing_taxon_states	; $href_hoh_info_of_infile_of_type->{supermatrix}{Nsequnconcat}{$taxon}++	; $href_hoh_info_of_infile_of_type->{$infile}{Ntaxa_mis}++	} #push @{$href_hol_seq_of_tax_of_infile->{$infile}}, ( $taxon, $string_missing_taxon_states	) }	# (2)
				elsif	( 'aa'		eq $href_hoh_info_of_infile_of_type->{$infile}{seqtype  }	){ $supermatrix_of_taxon{$taxon} .= $string_missing_taxon_states	; $href_hoh_info_of_infile_of_type->{supermatrix}{Nsequnconcat}{$taxon}++	; $href_hoh_info_of_infile_of_type->{$infile}{Ntaxa_mis}++	} #push @{$href_hol_seq_of_tax_of_infile->{$infile}}, ( $taxon, $string_missing_taxon_states	) }	# (3)
				else 	{ die "\n\tBUG-ERROR: Cannot concatenate sequence ", $taxon, " of file ", $infile, "!\n" }
			}
			############################################################ END CONCATENATION PROCESS
		}
		########################################################################################## END CONCATENATION OF EACH INFILE & SEQUENCE INFO EXTRACTION
		
		
		
		############################################################ START SUPERMATRIX TAXON ASSIGNMENT
		for my $taxon ( sort {$a<=>$b} keys %supermatrix_of_taxon ){
			
			push @{$href_hol_seq_of_tax_of_infile->{supermatrix}}, ( $taxon, $supermatrix_of_taxon{$taxon} ) ;
			
			# test print
			#if ( $taxon eq $href_hoh_info_of_infile_of_type->{supermatrix}{taxstruct} ){ print $taxon, $supermatrix_of_taxon{$taxon} } else{ print "Structur error!"}
			##
			
			$href_hoh_info_of_infile_of_type->{supermatrix}{Ntaxa}++ ;
			$href_hoh_info_of_infile_of_type->{supermatrix}{seqlength} = length $supermatrix_of_taxon{$taxon} ;
		}
		
		
		
		##############################
		## substitute supermatrix sequence end position by 1
		## supermatrix sequence endposition has been previously used as new startposition for next single infile
		$href_hoh_info_of_infile_of_type->{supermatrix}{seqend}-- ;
		##############################
		
		
		( %supermatrix_of_taxon ) = () ;
		############################################################ END SUPERMATRIX TAXON ASSIGNMENT
	}
	
	sub make_consensus{
		
		my $href_hol_seq_of_tax_of_infile	= $_[0] ;	# key: input file name; value ; list of taxa and corresponding sequences				-> IN (defined) / OUT (changed)
		my $sref_consensus_setup				= $_[1] ;	# defined consensus parameter (Freq, Maj, or Strict)									-> IN (defined) / OUT (unchanged)
		my $href_hoh_info_of_infile_of_type	= $_[2] ;	# key1: infilename, key2: type of info (e.g. 'seqlength'); info type associated value	-> IN (defined) / OUT (unchanged)
		my $href_name_list					= $_[3] ;	# key: taxon name; value : number of occurence in infiles								-> IN (defined) / OUT (changed)
		
		
		
		##############################
		# Definition of forward and reverse IUPAC ambiguity code
		my %iupac_code = (
						
						'forward' =>
						{
							'A'		=> 'A',
							'C'		=> 'C',
							'G'		=> 'G',
							'T'		=> 'T',
							'U'		=> 'T',
							'AC'	=> 'M',
							'AG'	=> 'R',
							'AT'	=> 'W',
							'CG'	=> 'S',
							'CT'	=> 'Y',
							'GT'	=> 'K',
							'ACG'	=> 'V',
							'ACT'	=> 'H',
							'AGT'	=> 'D',
							'CGT'	=> 'B',
							'ACGT'	=> 'N'
						},
						
						'reward' =>
						{
							'A'	=> 'A',
							'C'	=> 'C',
							'G'	=> 'G',
							'T'	=> 'T',
							'U'	=> 'T',
							'M'	=> 'AC',
							'R'	=> 'AG',
							'W'	=> 'AT',
							'S'	=> 'CG',
							'Y'	=> 'CT',
							'K'	=> 'GT',
							'V'	=> 'ACG',
							'H'	=> 'ACT',
							'D'	=> 'AGT',
							'B'	=> 'CGT',
							'N'	=> 'ACGT'
						},
		) ;
		##############################
		
		
		
		############################################################ START CONSENSUS BUILDING OF PREDEFINED SEQUENCE BLOCKS IN EACH FILE
		## For each infile...
		## ...identification of common sequence blocks which belong to the same taxon by analysing the taxon name prefix before the first underscore.
		## ...if taxon names consist of signs without underscore, the complete taxon name will be used to assign taxon name corresponding sequences
		## Sequences of the same defined block are stored in...
		## ...%hol_sequences_of_taxon-> key taxon name prefix with additional suffix ('_consensus'); value: list of sequences with identic taxon prefix
		## -----------------
		## e.g. File_1:
		## tax1_allel_0								->	ACGTTTTCGTTT...
		## tax1_allel_1								->	ACGTTTTAATTT...
		## tax1_allel_2								->	ACGTTTTGCTTT...
		## tax1										->	ACGTTTTTTTTT...
		## $consensus_seq_of_taxon{tax1_consensus}	= 	ACGTTTTNNTTT...
		## -----------------
		## Missing sequence information from taxon names with identic name prefix which are not included in the currently analysed infile are always considered during the consensus process !!
		## -----------------
		## e.g. File_1 without tax1_allel3:
		## tax1_allel_0								->	ACGTTTTCGTTT...
		## tax1_allel_1								->	ACGTTTTAATTT...
		## tax1_allel_2								->	ACGTTTTGCTTT...
		## tax1										->	ACGTTTTTTTTT...
		## tax1_allel_3								->	NNNNNNNNNNNN... -> missing sequence of tax1_allel3 in infile File_1 will be filled by N's during the concatenation process (if nucleotide data)
		## $consensus_seq_of_taxon{tax1_consensus}	= 	NNNNNNNNNNNN... -> strict consensus sequence (under option 'Strict')
		## -----------------
		%$href_name_list = () ;
		for my $infile ( sort {$a<=>$b} keys %$href_hol_seq_of_tax_of_infile ){
			
			
			##############################
			## Print file handling info
			print "\n\tGenerate ", $$sref_consensus_setup, " consensus blocks of ", $infile, "..." ;
			##############################
			
			
			
			##############################
			## Define hashvariabel in which consensus sequences of given infile are stored for each defined sequence block
			my %consensus_seq_of_taxon ;	# key: consensus taxon name; value: consensus sequence
			##############################
			
			
			
			##############################
			## Store identified sequences and associated taxon names of given infile in...
			## ...%sequence_of_taxon -> key: original taxon name; value: associated sequence
			my $aref_list_taxon_sequence_of_file	= exists( $href_hol_seq_of_tax_of_infile->{$infile}) ? \@{$href_hol_seq_of_tax_of_infile->{$infile}} :( ) ;
			my %sequence_of_taxon						= @$aref_list_taxon_sequence_of_file ;
			##############################
			
			
			
			############################################################ START IDENTIFICATION & ASSIGNMENT OF IDENTIFIED SEQUENCE BLOCKS
			## Identify number of taxa with identic name prefix (defined as one sequence block) and store each sequence as hashlist %hol_sequences_of_taxon
			## Ignore possible structure sequences and store them directly as consensus sequence in %consensus_seq_of_taxon
			my %hol_sequences_of_taxon ;
			for my $taxon ( sort keys %sequence_of_taxon ){
				
				
				##############################
				## Assign sequence blocks of same taxa of given input file and count assigned sequences
				## Store assigned sequences which have been identified by equal taxon prefix as hashlist in...
				## ...%hol_sequences_of_taxon -> key new taxon name with '_consensus' suffix; value: list of corresponding sequences
				## Possible structure sequences are excluded before consensus process
				unless ( $taxon eq $href_hoh_info_of_infile_of_type->{$infile}{taxstruct} ){
					
					my @name_parts	= split "_", $taxon ;
					my $con_taxon	= $name_parts[0]."_consensus" ;
					
					push @{$hol_sequences_of_taxon{$con_taxon}}, $sequence_of_taxon{$taxon} ;
				}
				##############################
				
				
				
				##############################
				## If structure sequence identified, push structure sequence in hashvariabel which stores each consensus sequence built for given infile
				else{ $consensus_seq_of_taxon{$taxon} = $sequence_of_taxon{$taxon} }
				##############################
			}
			############################################################ END IDENTIFICATION & ASSIGNMENT OF IDENTIFIED SEQUENCE BLOCKS
			
			
			
			############################################################ START CONSENSUS BUILDING OF DEFINED SEQUENCE BLOCKS OF GIVEN INFILE
			## Built consensus sequence of sequence pools with more than one defined sequence (in total overall infiles)
			##
			## Three different options to built a consensus sequence:
			##
			## (1) Most Frequent ('Freq') -> NUC & AA sequence data
			## ---------------------------------------------------------
			## Built consensus sequence by taking most frequent character state of each site
			## if two or more character states are equal frequent for a specific site take the iupac ambiguity code of these sites (nuc data)
			## or take '?' (amino acid data)
			## -----------------
			## e.g. File_1 (NUC):
			## tax1_allel_0								->	ACGTTTTCGTTT...
			## tax1_allel_1								->	ACGTTTTGTTTT...
			## tax1_allel_2								->	ACGTTTTGGTTT...
			## tax1										->	ACGTTTTTTTTT...
			## $consensus_seq_of_taxon{tax1_consensus}	= 	ACGTTTTGKTTT...
			##
			## e.g. File_1 (AA):
			## tax1_allel_0								->	ACGTTTTCGTTT...
			## tax1_allel_1								->	ACGTTTTATTTT...
			## tax1_allel_2								->	ACGTTTTGGTTT...
			## tax1										->	ACGTTTTTTTTT...
			## $consensus_seq_of_taxon{tax1_consensus}	= 	ACGTTTTG?TTT...
			## -----------------
			##
			## (2) Majority Frequent ('Maj') -> NUC & AA sequence data
			## ---------------------------------------------------------
			## Take character states which occure in more than 50% of sequences on a specific site
			## Otherwise take '?' for sites
			## -----------------
			## e.g. File_1 (NUC or AA):
			## tax1_allel_0								->	ACGTTTTGGTTT...
			## tax1_allel_1								->	ACGTTTTGTTTT...
			## tax1_allel_2								->	ACGTTTTGGTTT...
			## tax1										->	ACGTTTTTTTTT...
			## $consensus_seq_of_taxon{tax1_consensus}	= 	ACGTTTTG?TTT...
			## -----------------
			##
			## (3) Strict Consensus ('Strict') -> Only NUC (if defined for AA data skip to 'Freq')
			## ---------------------------------------------------------
			## Built a strict consensus sequence using iupac ambiguity code
			## -----------------
			## e.g. File_1 (NUC or AA):
			## tax1_allel_0								->	ACGTTTTGGTTT...
			## tax1_allel_1								->	ACGTTTTGTTTT...
			## tax1_allel_2								->	ACGTTTTGGTTT...
			## tax1										->	ACGTTTTTTTTT...
			## $consensus_seq_of_taxon{tax1_consensus}	= 	ACGTTTTKKTTT...
			## -----------------
			for my $taxon ( sort keys %hol_sequences_of_taxon ){
				
				
				##############################
				## Store sequences of sequence block in @$aref_sequence_pool
				## Count number of sequences states of sequence block
				## Count number of sequences of sequence block 
				my $aref_sequence_pool		= exists($hol_sequences_of_taxon{$taxon}) ? \@{$hol_sequences_of_taxon{$taxon}} :( ) ;
				my $N_states				= length $aref_sequence_pool->[0] ;
				my $N_pool_sequences		= @$aref_sequence_pool ;
				##############################
				
				
				
				############################## START GENERATION OF CONSENSUS SEQUENZ WITH DEFINED SEQUENCE BLOCKS > 1 SEQUENCE
				## Built consensus sequence from defined sequence pool (sequence block) if more than one sequence is defined in sequence pool of consensus taxon
				## Use substring to analyse each site position of a sequence one by one
				## Built consensus sequence stepwise until each sequence position has been analysed
				unless ( $N_pool_sequences == 1 ){
					
					
					##############################
					## define consensus sequence variable
					## which will be stepwise increased untile each sequence position has been analysed
					my $consensus_sequence	= () ;
					##############################
					
					
					
					############################################################ START ANALYSIS OF SINGLE SITE POSITIONS
					## Cout number of site states for each defined sequence of given sequence pool of given infile
					## Built consensus site state using either defined consensus option 'STRICT', 'Maj', or 'Freq'
					for my $seq_position ( 0 .. $N_states-1 ){
						
						
						############################################################ START IDENTIFICATION OF SITE STATE
						## Count character states at given sequence position
						my %counter_of_state ;
						for my $sequence_of_pool ( @$aref_sequence_pool ){
							
							my $character_state = substr ( $sequence_of_pool, $seq_position, 1 ) ;
							
							$counter_of_state{$character_state}++
						}
						##############################
						
						
						
						############################################################ START GENERATION OF SITE CONSENSUS 'Freq'
						## Built site specific consensus if 'Freq' option has been defined
						if 		(	$$sref_consensus_setup eq 'Freq' ){
							
							
							##############################
							## Create site consensus by taking the most frequent character state
							## If two or more character states are equally frequent...
							## ...take the corresponding iupac ambiguity if nucleotide data
							## ...take '?' if amino acid data
							my @sorted_occurence_of_states = sort{ $counter_of_state{$b}<=>$counter_of_state{$a} } keys %counter_of_state ;
							
							unless	( $sorted_occurence_of_states[1]																		){ $consensus_sequence .= $sorted_occurence_of_states[0] }
							elsif 	( $counter_of_state{$sorted_occurence_of_states[0]} > $counter_of_state{$sorted_occurence_of_states[1]}	){ $consensus_sequence .= $sorted_occurence_of_states[0] }
							elsif 	( $href_hoh_info_of_infile_of_type->{$infile}{seqtype} eq 'aa'											){ $consensus_sequence .= '?' }
							elsif	( $href_hoh_info_of_infile_of_type->{$infile}{seqtype} eq 'nu'											){
								
								
								##############################
								## (1) ...identify the equally highest frequent states and generate iupac ambiguity state
								## (2) ...states coded for missing data (?) or indel events (-) are not considered for ambiguity generation
								## (3) ...translate possible iupac ambiguity states to coding nucleotid estates & count nucleotide states
								## (4) ...sort list of counted nucleotide states
								## (5) ...transfer sorted list of nucleotdie states to string variabel
								## (6) ...translate string of nucleotides to corresponding iupac ambiguity code and enlarge consensus sequence
								my $first_most_freq_state = $sorted_occurence_of_states[0] ; my %counter_most_frequent_states ;
								my $gap = '-' ;
								
								for my $state ( 1 .. @sorted_occurence_of_states-1 ){
									
									# test print
									#print "\n", $sorted_occurence_of_states[$state], "\t", $counter_of_state{$sorted_occurence_of_states[$state]}, "\t", $counter_of_state{$first_most_freq_state},"\n";
									##
									
									if ( $counter_of_state{$sorted_occurence_of_states[$state]} == $counter_of_state{$first_most_freq_state} 	){ # (1)
										
										unless ( ( $state eq '?' ) || ( $state eq $gap ) ){	# (2)
											
											# (3)
											my		@states = split "", $iupac_code{reward}{$sorted_occurence_of_states[$state]} ;
											for (	@states ){ $counter_most_frequent_states{$_}++ }
										}
									}
								}
								
								# (3)
								my @states = split "", $iupac_code{reward}{$sorted_occurence_of_states[0]} ;
								for (	@states ){ $counter_most_frequent_states{$_}++ }
								####
								
								@sorted_occurence_of_states			 = sort keys %counter_most_frequent_states ;					# (4)
								my $string_most_frequent_states		 = join "", @sorted_occurence_of_states ;						# (5)
								
								if ( $iupac_code{forward}{$string_most_frequent_states} ){ $consensus_sequence .= $iupac_code{forward}{$string_most_frequent_states} }	# (6)
								else { $consensus_sequence .= '?' }
								
								# test print
								#print "\n", $infile, "\t", $taxon, " ", $consensus_sequence, " ", $string_most_frequent_states, "\n"; #exit;
								##
								
								( @sorted_occurence_of_states, %counter_most_frequent_states ) = () ;
								##############################
							}
							
							#print "\n", $infile, "\t", $href_hoh_info_of_infile_of_type->{$infile}{seqtype}, "\n" ;
							##############################
						}
						############################################################ END GENERATION OF SITE CONSENSUS 'Freq'
						
						
						
						############################################################ START GENERATION OF SITE CONSENSUS 'Maj'
						## Built site specific consensus if 'Maj' option has been defined
						elsif 	( $$sref_consensus_setup eq 'Maj' ){
							
							
							##############################
							## Create site consensus by taking the most frequent character state if number of characters
							## above 50%. Otherwise, take '?'
							## Determination by comparing number of most frequent character states with the sum of
							## remaining character states
							## Identic for nuc and aa
							my @sorted_occurence_of_states = sort{ $counter_of_state{$b}<=>$counter_of_state{$a} } keys %counter_of_state ;
							
							my $counter_remaining_states = 0 ;
							for ( 1 .. @sorted_occurence_of_states-1 ){ $counter_remaining_states += $counter_of_state{$sorted_occurence_of_states[$_]} }
							
							## Test Print
							#print "\n", $taxon, "\t", $seq_position, "\ts ", $counter_of_state{$sorted_occurence_of_states[0]}, "\th ", $counter_remaining_states, "\n";
							#####
							
							if ( $counter_of_state{$sorted_occurence_of_states[0]} > $counter_remaining_states )	{ $consensus_sequence .= $sorted_occurence_of_states[0] }
							else 																							{ $consensus_sequence .= '?' }
						}
						############################################################ END GENERATION OF SITE CONSENSUS 'Maj'
						
						
						
						############################################################ START GENERATION OF SITE CONSENSUS 'Strict'
						## Built site specific consensus if 'Strict' option has been defined
						## Create strict consensus of character states among sites by using iupac ambiguity code -> nuc
						## Create strict consensus of character states among sites by using 'X' as ambiguity state -> aa
						## sites with only missig and/or indel, and 1 kind of sequence state are summarized by the single character state
						## missing character states and indel states are kept unrecognized
						## site states which contain only indel and missing character states are summarized as '?'
						elsif 	( $$sref_consensus_setup eq 'Strict' ){
							
							
							##############################
							## Count the number of uninformative character states
							my @uninf_states	= ( '-', '?', 'X' ) ;
							my $counter_uninf	= 0 ;
							for ( @uninf_states ){ $counter_uninf += $counter_of_state{$_} }
							##############################
							
							
							
							##############################
							## Store informative character states as listvariabel
							my @states_found = keys %counter_of_state ;
							##############################
							
							
							
							##############################
							## If all pool sequences contain the same character state, take that character state
							if 		( $counter_of_state{$states_found[0]}	== $N_pool_sequences ){ $consensus_sequence .= $states_found[0] }
							##############################
							
							
							
							##############################
							## If total number of uninformative characters is equal the number of pool sequences take '?'
							elsif 	( $counter_uninf		== $N_pool_sequences ){ $consensus_sequence .= '?' }
							##############################
							
							
							
							##############################
							## If N of different nucleotide states > 1 -> Strict consensus of different informative site states
							## using IUPAC ambiguity code
							elsif ( $href_hoh_info_of_infile_of_type->{$infile}{seqtype} eq 'nu') {
								
								
								##############################
								## check for informative states
								## if site state of sequence pool is an ambiguity character
								## (1) split and count corresponding nucletotide states (e.g. $iupac_code{reward}{N} = ACGT -> @iupac_states( A C G T)
								## (2) sort list of counted states and take recoded ambiguity state (e.g. $iupac_code{forward}{ACGT} = N) as consensus state of that site
								%counter_of_state = () ;
								
								#### (1)
								for my $found_state ( @states_found ){
									
									if ( $iupac_code{reward}{$found_state} ){
										
										my @iupac_states = split "", $iupac_code{reward}{$found_state} ;
										for my $iupac_state ( @iupac_states ){ $counter_of_state{$iupac_state}++ }
									}
								}
								
								#### (2)
								my @recoded_states_found	= sort keys %counter_of_state ;
								my $recoded_state_string	= join "", @recoded_states_found ;
								
								## Test Print
								#print "\n", $taxon, "\t", $seq_position, "\ts ", $recoded_state_string, "\ti ", $iupac_code{forward}{$recoded_state_string}, "\n";
								#####
								
								$consensus_sequence .= $iupac_code{forward}{$recoded_state_string}
								##############################
							}
							##############################
							
							
							
							##############################
							## If N of different amino acid states > 1 -> Strict consensus of that site position = 'X'
							elsif ( $href_hoh_info_of_infile_of_type->{$infile}{seqtype} eq 'aa') { $consensus_sequence .= 'X' }
							##############################
							
							
							
							##############################
							## else die with a bug error prompt
							else{ die "\n\tBUG-ERROR: Cannot assign consensus state", $$sref_consensus_setup, "!\n" }
							##############################
						}
						############################################################ END GENERATION OF SITE CONSENSUS 'Strict'
						
						
						
						##############################
						## else die with a bug error prompt
						else{ die "\n\tBUG-ERROR: Cannot assign defined consensus parameter", $$sref_consensus_setup, "!\n" }
						##############################
						
						############################################################ END IDENTIFICATION OF SITE STATE
					}
					############################################################ END ANALYSIS OF SINGLE SITE POSITIONS
					
					#exit;
					
					##############################
					## Assign consensus sequence
					$consensus_seq_of_taxon{$taxon} = $consensus_sequence ;
					##############################
				}
				############################## END GENERATION OF CONSENSUS SEQUENZ WITH DEFINED SEQUENCE BLOCKS > 1 SEQUENCE
				
				
				
				############################## START TAKING SINGLE SEQUENCE IN SINGLE POOL AS CONSENSUS SEQUENCE
				else { $consensus_seq_of_taxon{$taxon} = $aref_sequence_pool->[0] }
				############################## END TAKING SINGLE SEQUENCE IN SINGLE POOL AS CONSENSUS SEQUENCE
			}
			############################################################ END CONSENSUS BUILDING OF PREDEFINED SEQUENCE BLOCKS OF GIVEN FILE
			
			
			
			##############################
			## Delete previous entry of taxon and sequences in given file
			## Clear taxon counter of infile due to compressed number of consensus sequences
			## Store consensus taxa name and consensus sequences as hashlist for given file
			## Count N consensus sequences of infile
			## Store new sequence names of all infiles as hashkey in %$href_name_list
			##############################
			$href_hol_seq_of_tax_of_infile	->{$infile} 		= () ;
			$href_hoh_info_of_infile_of_type->{$infile}{Ntaxa}	= 0 ;
			
			for my $con_taxon ( sort keys %consensus_seq_of_taxon ){
				
				push @{$href_hol_seq_of_tax_of_infile->{$infile}}, ( $con_taxon, $consensus_seq_of_taxon{$con_taxon} ) ;
				$href_hoh_info_of_infile_of_type->{$infile}{Ntaxa}++ ;
				
				$href_name_list->{$con_taxon}++
			}
			
			%consensus_seq_of_taxon = () ;
			##############################
		}
		############################################################ END CONSENSUS BUILDING OF PREDEFINED SEQUENCE BLOCKS IN EACH FILE
	}
	
	sub get_info{
		
		my $href_hol_seq_of_tax_of_infile		= $_[0] ;	# key: input file name; value ; list of taxa and corresponding sequences				-> IN (defined) / OUT (unchanged)
		my $href_hoh_info_of_infile_of_type		= $_[1] ;	# key1: infilename, key2: type of info (e.g. 'seqlength'); info type associated value	-> IN (defined) / OUT (changed)
		
		
		
		############################################################ START INFO STATE EXTRACTION OF CHARACTER STATES PER FILE
		for my $infile ( sort {$a<=>$b} keys %$href_hol_seq_of_tax_of_infile ){
			
			
			##############################
			## Store identified sequences and associated taxon names of given infile in...
			## ...%sequence_of_taxon -> key: original taxon name; value: associated sequence
			my $aref_list_taxon_sequence_of_file	= exists( $href_hol_seq_of_tax_of_infile->{$infile}) ? \@{$href_hol_seq_of_tax_of_infile->{$infile}} :( ) ;
			my %sequence_of_taxon						= @$aref_list_taxon_sequence_of_file ;
			##############################
			
			
			
			############################## START COUNTING SPECIFIC SITE STATES OF SEQUENCES IN INFILE
			## for each infile store...
			## (1) ...caharacter states of each site position in %informative_states_of_site_position
			my %states_of_site_position ; # key: site position number; value: list of character states found at given site
			## for each infile count...
			## (2) ...number of missing state characters -> nuc: '?'; aa: '?', 'X'
			## (3) ...number of indel states ('-')
			## (4) ...number of informative states -> nuc: 'A,C,G,T,U' aa: amino acid states
			## (4a)...number of GC states (only nucleotide data)
			## (5) ...number of ambiguity states -> only nuc: iupac ambiguity states
			for my $taxon ( sort keys %sequence_of_taxon ){
				
				
				############################################################ START CHECK OF DEFINED INFILE DATA OF GIVEN TAXON
				unless ( $taxon eq $href_hoh_info_of_infile_of_type->{$infile}{taxstruct} ){
					
					
					my @site_states			= split "", $sequence_of_taxon{$taxon} ;
					my $site_pos_counter	= 0 ;
					
					for my $state ( @site_states ){
						
						
						##############################
						## Count total number of infile states
						$href_hoh_info_of_infile_of_type->{$infile}{Nstates}++ ;
						##############################
						
						
						
						##############################
						## Store character states of given site position
						## as list in %informative_states_of_site_position
						push @{$states_of_site_position{$site_pos_counter}}, $state ;
						$site_pos_counter++;
						##############################
						
						
						
						##############################
						## If nucleodtide data
						if ( $href_hoh_info_of_infile_of_type->{$infile}{seqtype} eq 'nu' ){
							
							my $gap = '-' ;
							
							if 		( $state =~ /\?/						){ $href_hoh_info_of_infile_of_type->{$infile}{Nmissingstate	}++ }	# (2)
							elsif 	( $state =~ /$gap/						){ $href_hoh_info_of_infile_of_type->{$infile}{Ngapstate		}++ }	# (3)
							elsif	( $state =~ /A|C|G|T|U/					){
								
								$href_hoh_info_of_infile_of_type->{$infile}{$state}++ ; 
								$href_hoh_info_of_infile_of_type->{$infile}{infstate}++ ; 	# (4)
								
								if ( $state =~ /C|G/ ){ $href_hoh_info_of_infile_of_type->{$infile}{GC}++ }	# (4a)
							}
							elsif	( $state =~ /N|Y|R|W|S|K|M|D|V|H|B/		){ $href_hoh_info_of_infile_of_type->{$infile}{$state			}++; $href_hoh_info_of_infile_of_type->{$infile}{Nambstate	}++ }	# (5)
							else 	{ die "\n\tBUG-ERROR: Cannot assign nucleotide site state ", $state, " of taxon ", $taxon, " in file ", $infile, "!\n" }
						}
						##############################
						
						
						
						##############################
						## if amino acid data
						elsif ( $href_hoh_info_of_infile_of_type->{$infile}{seqtype} eq 'aa' ){
							
							my $gap = '-' ;
							
							if 		( $state =~ /\?|X/											){ $href_hoh_info_of_infile_of_type->{$infile}{Nmissingstate	}++ }	# (2)
							elsif 	( $state =~ /$gap/											){ $href_hoh_info_of_infile_of_type->{$infile}{Ngapstate		}++ }	# (3)
							elsif	( $state =~ /A|C|G|T|U|N|Y|R|W|S|K|M|D|V|H|B|F|L|I|P|Q|E/	){ $href_hoh_info_of_infile_of_type->{$infile}{$state			}++ ; $href_hoh_info_of_infile_of_type->{$infile}{infstate}++ }	# (4)
							else 	{ die "\n\tBUG-ERROR: Cannot assign amino acid site state ", $state, " of taxon ", $taxon, " in file ", $infile, "!\n" }
						}
						##############################
					}
				}
				############################################################ END CHECK OF DEFINED INFILE DATA OF GIVEN TAXON
			}
			############################## END COUNTING SPECIFIC SITE STATES OF SEQUENCES IN INFILE
			
			
			
			############################## START IDENTIFICATION OF PARSIMONY INFORMATIVE SITES
			## For each site position count number of informative site states stored in @{$states_of_site_position{$site_pos}}
			## If number of different informative site states > 2...
			## (1) ...set counter of informative site positions of given infile +1
			## (2) ...store informative site position number as list in @{$href_hoh_info_of_infile_of_type->{$infile}{list_inf_pos}}
			my $gap = '-' ;
			for my $site_pos ( sort {$a<=>$b} keys %states_of_site_position ){
				
				my %counter_inf_states ;
				my $aref_states = exists( $states_of_site_position{$site_pos} ) ? \@{$states_of_site_position{$site_pos}} :( ) ;
				
				for ( 0 .. @$aref_states-1 ){ 
					
					if ( $href_hoh_info_of_infile_of_type->{$infile}{seqtype} eq 'nu' ){
						
						if ( $aref_states->[$_] =~ /A|C|G|T|U/ ){ $counter_inf_states{$aref_states->[$_]}++ }
					}
					
					else{ unless ( $aref_states->[$_] =~ /\?|X|$gap/ ){ $counter_inf_states{$aref_states->[$_]}++ } }
				}
				
				 
				my $N_inf_states	= keys %counter_inf_states ;
				
				##############################
				## A site is parsimony-informative if it contains at least two types of nucleotides (or amino acids), 
				## and at least two of them occur with a minimum frequency of two
				if ( $N_inf_states >= 2 ){
					
					my @sort_states = sort { $counter_inf_states{$b} <=> $counter_inf_states{$a} } keys %counter_inf_states ;
					if ( ( $counter_inf_states{$sort_states[0]} > 1 ) && ( $counter_inf_states{$sort_states[1]} > 1 ) ){
						
								$href_hoh_info_of_infile_of_type->{$infile}{Ninf_sites	}++ ;			# (1)
						push @{	$href_hoh_info_of_infile_of_type->{$infile}{list_inf_pos}}, $site_pos	# (2)
					}
				}
				##############################
			}
			
			%states_of_site_position = () ;
			############################## END IDENTIFICATION OF PARSIMONY INFORMATIVE SITES
		}
		############################################################ END INFO STATE EXTRACTION OF CHARACTER STATES PER FILE
	}
	
	sub structure_handling{
		
		
		my $href_hol_seq_of_tax_of_infile	= $_[0] ;	# key: input file name; value ; list of taxa and corresponding sequences				-> IN (defined) / OUT (unchanged)
		my $href_hoh_info_of_infile_of_type	= $_[1] ;	# key1: infilename, key2: type of info (e.g. 'seqlength'); info type associated value	-> IN (defined) / OUT (changed)
		
		
		
		############################################################ START STRUCTURE ANALYSIS FOR EACH INFILE
		for my $infile ( sort {$a<=>$b} keys %$href_hol_seq_of_tax_of_infile ){
			
			
			if ( $href_hoh_info_of_infile_of_type->{$infile}{taxstruct} ){
				
				##############################
				## Store identified sequences and associated taxon names of given infile in...
				## ...%sequence_of_taxon -> key: original taxon name; value: associated sequence
				my $aref_list_taxon_sequence_of_file	= exists( $href_hol_seq_of_tax_of_infile->{$infile}) ? \@{$href_hol_seq_of_tax_of_infile->{$infile}} :( ) ;
				my %sequence_of_taxon						= @$aref_list_taxon_sequence_of_file ;
				##############################
				
				
				
				############################################################ START ASSIGNMENT OF LOOP & STEM REGIONS
				## Split taxon structure of given infile on each position and store single site structure states in @structures
				my @structures =  split ( "", $sequence_of_taxon{$href_hoh_info_of_infile_of_type->{$infile}{taxstruct}} ) ;
				#print "\n", $infile, "\t", @structures; 
				##############################
				
				
				
				##############################
				## Assign stem & loop regions by using site counter $i...
				## (1) ...if $i = '(', push @forward, $i
				## (2) ...if $i = ')', $j = pop @forward && push @pairs, ($j, $i)	-> collects paired stem positions
				## (3) ...if $i = '.', push @loops, $i								-> collects loop regions
				my  $i = 0 ;
				
				my (
						
						@forward,		# stores site numbers with open stem positions '('
						@loops,			# stores site numbers with closed stem positions ')'
						@pairs,			# stores pairing numbers open & closed stem positions joined by ':'
						@pairs_list,	# stores stem positions
				) ;
				
				for ( @structures ){ $i++ ;
					
					if 		( $_  =~ /\(/ ){								push @forward	, $i															}	# (1)
					elsif 	( $_  =~ /\)/ ){ my $pair_1 = pop @forward ;	push @pairs	, ( $pair_1.":".$i )	; push @pairs_list,( $pair_1, $i )	}	# (2)
					elsif 	( $_  =~ /\./ ){								push @loops	, $i															}	# (3)
					else 	{ die "\n\tBUG-ERROR: Cannot assign structer character ", $_, " in file ", $infile, "!\n"  }
				}
				
				#@pair_infos  =  reverse @pair_infos ;
				@pairs_list = sort {$a<=>$b} @pairs_list ;
				##############################
				
				############################################################ END ASSIGNMENT OF LOOP & STEM REGIONS
				
				
				
				############################################################ START EXTRACTION OF SECONDARY STRUCTURE INFO
				
				##############################
				## Store structure info...
				## (1) ...Number of structure positions
				## (2) ...Number of loop positions
				## (3) ...Number of stem positions
				## (4) ...Proportion of loop positions
				## (5) ...Proportion of stem positions
				## (6) ...String of stem positions e.g.: 1,4,6,7...
				## (7) ...String of loop positions e.g.: 2 3 5 8...
				$href_hoh_info_of_infile_of_type->{$infile}{Nstrucpos	}	= @structures ; #print "\n", $href_hoh_info_of_infile_of_type->{$infile}{Nstrucpos}; #exit;
				$href_hoh_info_of_infile_of_type->{$infile}{Nloops	}	= @loops ;
				$href_hoh_info_of_infile_of_type->{$infile}{Nstems	}	= $href_hoh_info_of_infile_of_type->{$infile}{Nstrucpos} - $href_hoh_info_of_infile_of_type->{$infile}{Nloops} ;
				$href_hoh_info_of_infile_of_type->{$infile}{Ploops	}	= sprintf "%.3f", ( ( $href_hoh_info_of_infile_of_type->{$infile}{Nloops} / $href_hoh_info_of_infile_of_type->{$infile}{Nstrucpos} ) * 100 );
				$href_hoh_info_of_infile_of_type->{$infile}{Pstems	}	= sprintf "%.3f", ( $href_hoh_info_of_infile_of_type->{$infile}{Nstems} / $href_hoh_info_of_infile_of_type->{$infile}{Nstrucpos} ) * 100 ;
				$href_hoh_info_of_infile_of_type->{$infile}{Stems		}	= join " ", @pairs ;
				$href_hoh_info_of_infile_of_type->{$infile}{Pairs		}	= join ",", @pairs ;
				$href_hoh_info_of_infile_of_type->{$infile}{Pairlist	}	= join ",", @pairs_list ;
				$href_hoh_info_of_infile_of_type->{$infile}{Loops		}	= join " ", @loops ;
				##############################
				
				############################################################ END EXTRACTION OF SECONDARY STRUCTURE INFO
			}
		}
		############################################################ END STRUCTURE ANALYSIS FOR EACH INFILE
	}
	
	sub print_out{
		
		
		my $href_hol_seq_of_tax_of_infile		= $_[0] ;	# key: input file name; value ; list of taxa and corresponding sequences				-> IN (defined) / OUT (unchanged)
		my $href_hoh_info_of_infile_of_type		= $_[1] ;	# key1: infilename, key2: type of info (e.g. 'seqlength'); info type associated value	-> IN (defined) / OUT (unchanged)
		my $href_outfile_name					= $_[2] ;	# Outfilename of output option															-> IN (defined) / OUT (unchanged)
		my $aref_parameter_all					= $_[3] ;	# List of parameter options of data concatenation										-> IN (defined) / OUT (unchanged)
		my $aref_parameter_inf					= $_[4] ;	# List of parameter options of info print out											-> IN (defined) / OUT (unchanged)
		my $aref_parameter_phy					= $_[5] ;	# List of parameter options of PHYLIP print OUT											-> IN (defined) / OUT (unchanged)
		my $aref_parameter_nex					= $_[6] ;	# List of parameter options of NEXUS print OUT											-> IN (defined) / OUT (unchanged)
		my $aref_parameter_con					= $_[7] ;	# List of consensus options																-> IN (defined) / OUT (unchanged)
		my $aref_parameter_fil					= $_[8] ;	# List of of file handling																-> IN (defined) / OUT (unchanged)
		my $aref_parameter_fas					= $_[9] ;	# defined parameter option of FASTA print OUT											-> IN (defined) / OUT (unchanged)
		my $aref_parameter_3rd					= $_[10];	# List of third position handling														-> IN (defined) / OUT (unchanged)
		my $aref_parameter_ryc					= $_[11];	# List of RY coding																		-> IN (defined) / OUT (unchanged)
		my $aref_parameter_tra					= $_[12];	# List of sequence translation options													-> IN (defined) / OUT (unchanged)
		my $aref_parameter_par					= $_[13];	# Print parsimonious sites as extra msa file											-> IN (defined) / OUT (unchanged)
		my $aref_parameter_ren					= $_[14];	# Rename taxon names of given ifiles													-> IN (defined) / OUT (unchanged)
		my $aref_parameter_prt					= $_[15];	# Print partition files for concatenated data											-> IN (defined) / OUT (unchanged)
		my $aref_parameter_pro					= $_[16];	# Start prottest analyses for aa data													-> IN (defined) / OUT (unchanged)
		my $aref_parameter_mis					= $_[17];	# Replacement code of missing gene sequences 											-> IN (defined) / OUT (unchanged)
		my $sref_structure_seq					= $_[18];	# defined if structure sequence present in infiles										-> IN (defined) / OUT (unchanged)
		
		
		
		############################################################ START SUPERMATRIX PRINT OUT
		
		##############################
		## define output file names according to chosen infile handling option
		my @outfiles ;
		if 		( $aref_parameter_fil->[0] eq 'Supermatrix'			)	{ @outfiles = 'supermatrix' }
		else 															{ @outfiles = sort keys %$href_hol_seq_of_tax_of_infile }
		##############################
		
		
		
		##############################
		## OUTPUT -> FASTA format 
		unless ( $aref_parameter_fas->[0] eq 'NO' ){
			
			&fas_out(
							\@outfiles,							# list of filenames for print OUT														-> IN (defined) / OUT (unchanged)
							\%$href_hol_seq_of_tax_of_infile,	# key: input file name; value ; list of taxa and corresponding sequences				-> IN (defined) / OUT (unchanged)
							\%$href_hoh_info_of_infile_of_type,	# key1: infilename, key2: type of info (e.g. 'seqlength'); info type associated value	-> IN (defined) / OUT (changed)
							\'non_parsimony'					# if 'parsimony' do not set conversion counter +1										-> IN (defined) / OUT (unchanged)
			) ;
			
			if ( ( $aref_parameter_prt->[0] eq 'Supermatrix' ) && ( $aref_parameter_fil->[0] =~ /Supermatrix/ ) ){
				
				&raxml_out_prt(
							\$aref_parameter_prt->[0],			# Print partition files for concatenated data											-> IN (defined) / OUT (unchanged)
							\%$href_hol_seq_of_tax_of_infile,	# key: input file name; value ; list of taxa and corresponding sequences				-> IN (defined) / OUT (unchanged)
							\%$href_hoh_info_of_infile_of_type,	# key1: infilename, key2: type of info (e.g. 'seqlength'); info type associated value	-> IN (defined) / OUT (unchanged)
							\$aref_parameter_pro->[0],			# parameter for prottest analyses														-> IN (defined) / OUT (unchanged)
				) ;
			}
		}
		##############################
		
		
		
		##############################
		## Additional output (phylip STRICT or RELAXED)
		unless ( $aref_parameter_phy->[0] eq 'NO' ){
			
			&phylip_out(
							\@outfiles,							# list of filenames for print OUT														-> IN (defined) / OUT (unchanged)
							\%$href_hol_seq_of_tax_of_infile,	# key: input file name; value ; list of taxa and corresponding sequences				-> IN (defined) / OUT (unchanged)
							\$aref_parameter_phy->[0],		# phylip output parameter ('strict' or 'relaxed')										-> IN (defined) / OUT (unchanged)
							\%$href_hoh_info_of_infile_of_type,	# key1: infilename, key2: type of info (e.g. 'seqlength'); info type associated value	-> IN (defined) / OUT (changed)
							\'non_parsimony'					# if 'parsimony' do not set conversion counter +1										-> IN (defined) / OUT (unchanged)
			) ;
			
			if ( ( $aref_parameter_prt->[0] eq 'Supermatrix' ) && ( $aref_parameter_fil->[0] =~ /Supermatrix/ ) ){
				
				&raxml_out_prt(
							\$aref_parameter_prt->[0],			# Print partition files for concatenated data											-> IN (defined) / OUT (unchanged)
							\%$href_hol_seq_of_tax_of_infile,	# key: input file name; value ; list of taxa and corresponding sequences				-> IN (defined) / OUT (unchanged)
							\%$href_hoh_info_of_infile_of_type,	# key1: infilename, key2: type of info (e.g. 'seqlength'); info type associated value	-> IN (defined) / OUT (unchanged)
							\$aref_parameter_pro->[0],			# parameter for prottest analyses														-> IN (defined) / OUT (unchanged)
				) ;
			}
		}
		##############################
		
		
		
		##############################
		## Additional output (NEXUS 'BLOCK' or 'MrBAYES')
		unless ( $aref_parameter_nex->[0] eq 'NO' ){
			
			&nexus_out(
							\@outfiles,							# list of filenames for print OUT														-> IN (defined) / OUT (unchanged)
							\$aref_parameter_nex->[0],		# nexus output parameter ('BLOCK' or 'MrBAYES')											-> IN (defined) / OUT (unchanged)
							\%$href_hol_seq_of_tax_of_infile,	# key: input file name; value ; list of taxa and corresponding sequences				-> IN (defined) / OUT (unchanged)
							\%$href_hoh_info_of_infile_of_type,	# key1: infilename, key2: type of info (e.g. 'seqlength'); info type associated value	-> IN (defined) / OUT (unchanged)
							\'non_parsimony'					# if 'parsimony' do not set conversion counter +1										-> IN (defined) / OUT (unchanged)
			) ;
			
			if ( ( $aref_parameter_prt->[0] eq 'Supermatrix' ) && ( $aref_parameter_fil->[0] =~ /Supermatrix/ ) ){
				
				&nexus_out_prt(
							\$aref_parameter_nex->[0],		# nexus output parameter ('BLOCK' or 'MrBAYES')											-> IN (defined) / OUT (unchanged)
							\$aref_parameter_prt->[0],		# # Print partition files for concatenated data											-> IN (defined) / OUT (unchanged)
							\%$href_hol_seq_of_tax_of_infile,	# key: input file name; value ; list of taxa and corresponding sequences				-> IN (defined) / OUT (unchanged)
							\%$href_hoh_info_of_infile_of_type,	# key1: infilename, key2: type of info (e.g. 'seqlength'); info type associated value	-> IN (defined) / OUT (unchanged)
				) ;
			}
		}
		##############################
		
		
		
		##############################
		## Additional output (PARSIMONY INF SITES)
		unless ( $aref_parameter_par->[0] eq 'NO' ){
			
			&parsimony_print(
							\@outfiles,								# list of filenames for print OUT														-> IN (defined) / OUT (unchanged)
							\%$href_hol_seq_of_tax_of_infile,		# key: input file name; value ; list of taxa and corresponding sequences				-> IN (defined) / OUT (unchanged)
							\%$href_hoh_info_of_infile_of_type,		# key1: infilename, key2: type of info (e.g. 'seqlength'); info type associated value	-> IN (defined) / OUT (unchanged)
							\$aref_parameter_fas->[0],			# defined parameter option of FASTA print OUT											-> IN (defined) / OUT (unchanged)
							\$aref_parameter_nex->[0],			# defined parameter option of NEXUS print OUT											-> IN (defined) / OUT (unchanged)
							\$aref_parameter_phy->[0],			# defined parameter option of PHYLIP print OUT											-> IN (defined) / OUT (unchanged)
			) ;
		}
		############################################################ END SUPERMATRIX PRINT OUT
		
		
		
		############################################################ START INFO PRINT OUT
		
		##############################
		## Print basal info
		if ( $aref_parameter_inf->[0] eq 'NO' ){
			
			&info_small (
							\%$href_hol_seq_of_tax_of_infile,		# key: input file name; value ; list of taxa and corresponding sequences				-> IN (defined) / OUT (unchanged)
							\%$href_hoh_info_of_infile_of_type,		# key1: infilename, key2: type of info (e.g. 'seqlength'); info type associated value	-> IN (defined) / OUT (unchanged)
							\$href_outfile_name->{info},			# info outfile name																		-> IN (defined) / OUT (unchanged)
							\@$aref_parameter_all,					# defined parameter option of data concatenation										-> IN (defined) / OUT (unchanged)
							\@$aref_parameter_inf,					# defined parameter option of info print out											-> IN (defined) / OUT (unchanged)
							\@$aref_parameter_phy,					# defined parameter option of PHYLIP print OUT											-> IN (defined) / OUT (unchanged)
							\@$aref_parameter_nex,					# defined parameter option of NEXUS print OUT											-> IN (defined) / OUT (unchanged)
							\@$aref_parameter_con,					# defined consensus option																-> IN (defined) / OUT (unchanged)
							\@$aref_parameter_fil,					# defined parameter option of file handling												-> IN (defined) / OUT (unchanged)
							\@$aref_parameter_fas,					# defined parameter option of FASTA print OUT											-> IN (defined) / OUT (unchanged)
							\@$aref_parameter_3rd,					# List of third position handling														-> IN (defined) / OUT (unchanged)
							\@$aref_parameter_ryc,					# List of RY coding																		-> IN (defined) / OUT (unchanged)
							\@$aref_parameter_tra,					# List of Sequence transaltion options													-> IN (defined) / OUT (unchanged)
							\@$aref_parameter_par,					# Print parsimonious sites as extra msa file											-> IN (defined) / OUT (unchanged)
							\@$aref_parameter_ren,					# Rename taxon names of given ifiles													-> IN (defined) / OUT (unchanged)
							\@$aref_parameter_prt,					# Print partition files for concatenated data											-> IN (defined) / OUT (unchanged)
							\@$aref_parameter_pro,					# Start prottest analyses for aa data													-> IN (defined) / OUT (unchanged)
							\@$aref_parameter_mis,					# Replacement code of missing gene sequences 											-> IN (defined) / OUT (unchanged)
			)
		}
		##############################
		
		
		
		##############################
		## Print complete Info
		elsif ( $aref_parameter_inf->[0] eq 'YES' ){
			
			&info_all (
							\%$href_hol_seq_of_tax_of_infile,		# key: input file name; value ; list of taxa and corresponding sequences				-> IN (defined) / OUT (unchanged)
							\%$href_hoh_info_of_infile_of_type,		# key1: infilename, key2: type of info (e.g. 'seqlength'); info type associated value	-> IN (defined) / OUT (unchanged)
							\$href_outfile_name->{info},			# info outfile name																		-> IN (defined) / OUT (unchanged)
							\$href_outfile_name->{structure},		# structure sequence outfile name														-> IN (defined) / OUT (unchanged)
							\@$aref_parameter_all,					# defined parameter option of data concatenation										-> IN (defined) / OUT (unchanged)
							\@$aref_parameter_inf,					# defined parameter option of info print out											-> IN (defined) / OUT (unchanged)
							\@$aref_parameter_phy,					# defined parameter option of PHYLIP print OUT											-> IN (defined) / OUT (unchanged)
							\@$aref_parameter_nex,					# defined parameter option of NEXUS print OUT											-> IN (defined) / OUT (unchanged)
							\@$aref_parameter_con,					# defined consensus option																-> IN (defined) / OUT (unchanged)
							\@$aref_parameter_fil,					# defined parameter option of file handling												-> IN (defined) / OUT (unchanged)
							\@$aref_parameter_fas,					# defined parameter option of FASTA print OUT											-> IN (defined) / OUT (unchanged)
							\@$aref_parameter_3rd,					# List of third position handling														-> IN (defined) / OUT (unchanged)
							\@$aref_parameter_ryc,					# List of RY coding																		-> IN (defined) / OUT (unchanged)
							\@$aref_parameter_tra,					# List of Sequence transaltion options													-> IN (defined) / OUT (unchanged)
							\@$aref_parameter_par,					# Print parsimonious sites as extra msa file											-> IN (defined) / OUT (unchanged)
							\@$aref_parameter_ren,					# Rename taxon names of given ifiles													-> IN (defined) / OUT (unchanged)
							\@$aref_parameter_prt,					# Print partition files for concatenated data											-> IN (defined) / OUT (unchanged)
							\@$aref_parameter_pro,					# Start prottest analyses for aa data													-> IN (defined) / OUT (unchanged)
							\@$aref_parameter_mis,					# Replacement code of missing gene sequences 											-> IN (defined) / OUT (unchanged)
							\$$sref_structure_seq,					# defined if structure sequence present in infiles										-> IN (defined) / OUT (unchanged)
			)
		}
		##############################
		
		
		
		##############################
		## Print BUG RPORT
		else{ die "\n\tBUG-ERROR: Cannot assign info option ", $aref_parameter_inf->[0], "!\n\tPlease, report BUG to system developer!" }
		##############################
		
		
		
		########################################################################################## START SUBROUTINES OF &print_out
		
		sub info_small{
			
			my $href_hol_seq_of_tax_of_infile	= $_[0] ;	# key: input file name; value ; list of taxa and corresponding sequences				-> IN (defined) / OUT (unchanged)
			my $href_hoh_info_of_infile_of_type	= $_[1] ;	# key1: infilename, key2: type of info (e.g. 'seqlength'); info type associated value	-> IN (defined) / OUT (unchanged)
			my $sref_output_name				= $_[2] ;	# info outfile name																		-> IN (defined) / OUT (unchanged)
			my $aref_parameter_all				= $_[3] ;	# defined parameter option of data concatenation										-> IN (defined) / OUT (unchanged)
			my $aref_parameter_inf				= $_[4] ;	# defined parameter option of info print out											-> IN (defined) / OUT (unchanged)
			my $aref_parameter_phy				= $_[5] ;	# defined parameter option of PHYLIP print OUT											-> IN (defined) / OUT (unchanged)
			my $aref_parameter_nex				= $_[6] ;	# defined parameter option of NEXUS print OUT											-> IN (defined) / OUT (unchanged)
			my $aref_parameter_con				= $_[7] ;	# defined consensus option																-> IN (defined) / OUT (unchanged)
			my $aref_parameter_fil				= $_[8] ;	# defined parameter option of file handling												-> IN (defined) / OUT (unchanged)
			my $aref_parameter_fas				= $_[9] ;	# defined parameter option of FASTA print OUT											-> IN (defined) / OUT (unchanged)
			my $aref_parameter_3rd				= $_[10];	# List of third position handling														-> IN (defined) / OUT (unchanged)
			my $aref_parameter_ryc				= $_[11];	# List of RY coding																		-> IN (defined) / OUT (unchanged)
			my $aref_parameter_tra				= $_[12];	# List of sequence translation options													-> IN (defined) / OUT (unchanged)
			my $aref_parameter_par				= $_[13];	# Print parsimonious sites as extra msa file											-> IN (defined) / OUT (unchanged)
			my $aref_parameter_ren				= $_[14];	# Rename taxon names of given ifiles													-> IN (defined) / OUT (unchanged)
			my $aref_parameter_prt				= $_[15];	# Print partition files for concatenated data											-> IN (defined) / OUT (unchanged)
			my $aref_parameter_pro				= $_[16];	# Start prottest analyses for aa data													-> IN (defined) / OUT (unchanged)
			my $aref_parameter_mis				= $_[17];	# Replacement code of missing gene sequences 											-> IN (defined) / OUT (unchanged)
			
			
			##############################
			## Open OUT info file & Print info header
			open  OUT_info, ">$$sref_output_name" ||  warn "\n\t!FILE-ERROR!: Cannot OPEN OUT ", $$sref_output_name, "!\n" ;
			print OUT_info "FASconCAT INFO:";
			##############################
			
			
			
			##############################
			## If supermatrix is printed out
			## (1) ...print infile name and supermatrix range of each infile
			## (2) ...print out output format of supermatrix
			if ( $aref_parameter_fil->[0] =~ /Supermatrix/ ){
				
				# (1)
				print OUT_info "\nFile\tSTART POSITION SUPERMATRIX\tEND POSITION SUPERMATRIX\n" ;
				
				for my $infile ( sort keys %$href_hol_seq_of_tax_of_infile ){
					
					( my $file = $infile ) =~ s/.fas$|.aln$|.FASTA$|.phy$// ;
					
					print OUT_info	$file, "\t", $href_hoh_info_of_infile_of_type->{$infile}{seqstart}, "\t", $href_hoh_info_of_infile_of_type->{$infile}{seqend}, "\n"
				}
				
				
				# (2)
				print OUT_info	"\nSUPERMATRIX FILE\tPRINTED AS FASTA\tPRINTED AS NEXUS\tPRINTED AS PHYLIP\nsupermatrix";
				
				( $href_hoh_info_of_infile_of_type->{supermatrix}{convertedfas} == 1 ) ? ( print OUT_info "\tYES" ) : ( print OUT_info "\tNO" ) ;
				( $href_hoh_info_of_infile_of_type->{supermatrix}{convertednex} == 1 ) ? ( print OUT_info "\tYES" ) : ( print OUT_info "\tNO" ) ;
				( $href_hoh_info_of_infile_of_type->{supermatrix}{convertedphy} == 1 ) ? ( print OUT_info "\tYES" ) : ( print OUT_info "\tNO" ) ;
			}
			##############################
			
			
			
			##############################
			## If single data sets have been printed out in converted file format
			## print filename & output formats
			## Only for input files, not for supermatrix output
			if ( $aref_parameter_fil->[0] =~ /Convert/ ){
				
				print OUT_info "\nFile\tL SEQUENCE [bp]\tCONVERTED TO FASTA\tCONVERTED TO NEXUS\tCONVERTED TO PHYLIP" ;
				
				for my $infile ( sort {$a<=>$b} keys %$href_hol_seq_of_tax_of_infile ){
					
					unless ( $infile eq 'supermatrix'){
						
						( my $file = $infile ) =~ s/.fas$|.aln$|.FASTA$|.phy$// ;
						
						print OUT_info "\n", $file, "\t", $href_hoh_info_of_infile_of_type->{$infile}{seqlength} ;
						
						( $href_hoh_info_of_infile_of_type->{$infile}{convertedfas} == 1 ) ? ( print OUT_info "\tYES" ) : ( print OUT_info "\tNO" ) ;
						( $href_hoh_info_of_infile_of_type->{$infile}{convertednex} == 1 ) ? ( print OUT_info "\tYES" ) : ( print OUT_info "\tNO" ) ;
						( $href_hoh_info_of_infile_of_type->{$infile}{convertedphy} == 1 ) ? ( print OUT_info "\tYES" ) : ( print OUT_info "\tNO" ) ;
					}
				}
			}
			##############################
			
			
			
			##############################
			## Print... 
			## (1) ...N of taxon concatenations (total)
			## (2) ...N informative taxon sequences
			## (3) ...N missing taxon sequences
			if ( $aref_parameter_fil->[0] =~ /Supermatrix/ ){
				print OUT_info "\n",
								"\nTaxon check\n",
								"Taxon\tN CONCATENATED SEQs\tN INFORMATIVE SEQs\tN MISSING SEQs\n"
				;
				
				my %supermatrix_seq_of_taxon = @{$href_hol_seq_of_tax_of_infile->{supermatrix}} ;
				for my $taxon ( sort {$a<=>$b} keys %supermatrix_seq_of_taxon ){
					
					print OUT_info $taxon,
												"\t", $href_hoh_info_of_infile_of_type->{supermatrix}{taxconcat		}{$taxon},	# (1)
												"\t", $href_hoh_info_of_infile_of_type->{supermatrix}{Nseqconcat	}{$taxon},	# (2)
												"\t", $href_hoh_info_of_infile_of_type->{supermatrix}{Nsequnconcat	}{$taxon},	# (3)
												"\n"
				}
			}
			##############################
			
			
			
			##############################
			## Terminal Info Print...
			## (1) ...Site Ranges of given infiles in the concatenated supermatrix (if any...)
			## (2) ...the concatenated supermatrix (if given)
			
			# (1)
			if ( $aref_parameter_fil->[0] =~ /Supermatrix/ ){
				
				print "\n\n\tSupermatrix Site Range Of File...\n\t" ;
				for my $infile ( sort keys %$href_hol_seq_of_tax_of_infile ){
					
					( my $file = $infile ) =~ s/.fas$|.aln$|.FASTA$|.phy$// ;
					
					unless ( $file eq 'supermatrix' ){ print $file, ":\tSite Pos.", $href_hoh_info_of_infile_of_type->{$infile}{seqstart}, " to ", $href_hoh_info_of_infile_of_type->{$infile}{seqend}, "\n\t" }
				}
			}
			
			# (2)
			if ( $aref_parameter_fil->[0] =~ /Supermatrix/ ){
				
				if ( $href_hoh_info_of_infile_of_type->{supermatrix}{Nmissingstate} ){
					
					print "\n\tSupermatrix:\t", ( sprintf "%.3f", ( $href_hoh_info_of_infile_of_type->{supermatrix}{Nmissingstate} / $href_hoh_info_of_infile_of_type->{supermatrix}{Nstates} ) * 100 ), "%"
				}
				
				else{ print "\n\tSupermatrix:\t", 0, "%" }
			}
			##############################
			
			
			
			##############################
			# Print info header
			print OUT_info "\n"							,
							"\n------------------------------------------------------------" ,
							"\nDefined Options:"			,
							"\n-----------------"			,
							"\nREAD IN ALL files:\t"		,	$aref_parameter_all->[0],
							"\nREAD IN SINGLE files:\t"		,	$aref_parameter_all->[1],
							"\nFILE HANDLING:\t"			,	$aref_parameter_fil->[0],
							"\nSEQUENCE TRANSLATION\t"		,	$aref_parameter_tra->[0],
							"\nEXCLUDING 3rd POSITION\t"	,	$aref_parameter_3rd->[0],
							"\nRY CODING\t"					,	$aref_parameter_ryc->[0],
							"\nBUILD CONSENSUS\t"			,	$aref_parameter_con->[0],
							"\nRENAME SEQUENCE NAMES\t"		,	$aref_parameter_ren->[0],
							"\nABSENT SEQUENCE CODING\t"	,	$aref_parameter_mis->[0],
							"\n"							,
							"\n-----------------"			,
							"\nOUTPUT FORMAT:"				,
							"\nPHYLIP:\t"					,	$aref_parameter_phy->[0],
							"\nNEXUS:\t"					,	$aref_parameter_nex->[0],
							"\nFASTA:\t"					,	$aref_parameter_fas->[0],
							"\nPARTITION:\t"				,	$aref_parameter_prt->[0],
							"\nPARSIMONY:\t"				,	$aref_parameter_par->[0],
							"\nPROTTEST:\t"					,	$aref_parameter_pro->[0],
							"\n------------------------------------------------------------\n"
			;
			##############################
			
			
			
			##############################
			## close info file
			close OUT_info ;
			##############################
		}
		
		sub info_all{
			
			
			my $href_hol_seq_of_tax_of_infile		= $_[0] ;	# key: input file name; value ; list of taxa and corresponding sequences				-> IN (defined) / OUT (unchanged)
			my $href_hoh_info_of_infile_of_type		= $_[1] ;	# key1: infilename, key2: type of info (e.g. 'seqlength'); info type associated value	-> IN (defined) / OUT (unchanged)
			my $sref_output_name					= $_[2] ;	# info outfile name																		-> IN (defined) / OUT (unchanged)
			my $sref_output_name_structure			= $_[3] ;	# structure sequence outfile name														-> IN (defined) / OUT (unchanged)
			my $aref_parameter_all					= $_[4] ;	# defined parameter option of data concatenation										-> IN (defined) / OUT (unchanged)
			my $aref_parameter_inf					= $_[5] ;	# defined parameter option of info print out											-> IN (defined) / OUT (unchanged)
			my $aref_parameter_phy					= $_[6] ;	# defined parameter option of PHYLIP print OUT											-> IN (defined) / OUT (unchanged)
			my $aref_parameter_nex					= $_[7] ;	# defined parameter option of NEXUS print OUT											-> IN (defined) / OUT (unchanged)
			my $aref_parameter_con					= $_[8] ;	# defined consensus option																-> IN (defined) / OUT (unchanged)
			my $aref_parameter_fil					= $_[9] ;	# defined parameter option of file handling												-> IN (defined) / OUT (unchanged)
			my $aref_parameter_fas					= $_[10];	# defined parameter option of FASTA print OUT											-> IN (defined) / OUT (unchanged)
			my $aref_parameter_3rd					= $_[11];	# List of third position handling														-> IN (defined) / OUT (unchanged)
			my $aref_parameter_ryc					= $_[12];	# List of RY coding																		-> IN (defined) / OUT (unchanged)
			my $aref_parameter_tra					= $_[13];	# List of sequence translation options													-> IN (defined) / OUT (unchanged)
			my $aref_parameter_par					= $_[14];	# Print parsimonious sites as extra msa file											-> IN (defined) / OUT (unchanged)
			my $aref_parameter_ren					= $_[15];	# Rename taxon names of given ifiles													-> IN (defined) / OUT (unchanged)
			my $aref_parameter_prt					= $_[16];	# Print partition files for concatenated data											-> IN (defined) / OUT (unchanged)
			my $aref_parameter_pro					= $_[17];	# Start prottest analyses for aa data													-> IN (defined) / OUT (unchanged)
			my $aref_parameter_mis					= $_[18];	# Replacement code of missing gene sequences 											-> IN (defined) / OUT (unchanged)
			my $structure_seq						= $_[19];	# defined if structure sequence present in infiles										-> IN (defined) / OUT (unchanged)
			
			
			##############################
			## Open OUT info file
			open  OUT_info, ">$$sref_output_name" ||  warn "\n\t!FILE-ERROR!: Cannot OPEN OUT ", $$sref_output_name, "!\n" ;
			print OUT_info "FASconCAT INFO:\n";
			##############################
			
			
			
			############################################################ START SEQUENCE INFO PRINT OUT
			## If supermatrix is printed out
			## (1) ...print out output format of supermatrix and additional supermatrix info
			## (2) ...print infile name and supermatrix range of each infile and sequence infos
			if ( $aref_parameter_fil->[0] =~ /Supermatrix/ ){
				
				# (1)
				print OUT_info	"SUPERMATRIX FILE\tP CODING CHARACTERS [%]\tP MISSING DATA [?] [%]\tP INDELS [%]\tP AMBIGUITIES [%]\tP PARSIMONY INFORMATIVE SITES [%]\tP GC CONTENT [%]\tPRINTED AS FASTA\tPRINTED AS NEXUS\tPRINTED AS PHYLIP\nsupermatrix";
				
				my @printouts ;
				for my $state ( qw/infstate Nmissingstate Ngapstate Nambstate/ ){
					
					if ( $href_hoh_info_of_infile_of_type->{supermatrix}{$state} ){
						
						print OUT_info "\t", ( sprintf "%.3f", ( $href_hoh_info_of_infile_of_type->{supermatrix}{$state} / $href_hoh_info_of_infile_of_type->{supermatrix}{Nstates} ) * 100 )
					}
					
					else{ print OUT_info "\t", 0 }
				}
				
				if ( $href_hoh_info_of_infile_of_type->{supermatrix}{infstate} ){
					
					print OUT_info "\t", ( sprintf "%.3f", ( $href_hoh_info_of_infile_of_type->{supermatrix}{Ninf_sites} / $href_hoh_info_of_infile_of_type->{supermatrix}{seqlength} ) * 100 )
				}
				else{ print OUT_info "\t", 0 }
				
				if ( $href_hoh_info_of_infile_of_type->{supermatrix}{seqtype} eq 'nu' ){
					
					if ( $href_hoh_info_of_infile_of_type->{supermatrix}{GC} ){
						
						print OUT_info "\t", ( sprintf "%.3f", ( $href_hoh_info_of_infile_of_type->{supermatrix}{GC} / $href_hoh_info_of_infile_of_type->{supermatrix}{Nstates} ) * 100 )
					}
					
					else{ print OUT_info "\t", 0 }
				}
				else{  print OUT_info "\tonly nuc data" }
				
				( $href_hoh_info_of_infile_of_type->{supermatrix}{convertedfas} == 1 ) ? ( print OUT_info "\tYES" ) : ( print OUT_info "\tNO" ) ;
				( $href_hoh_info_of_infile_of_type->{supermatrix}{convertednex} == 1 ) ? ( print OUT_info "\tYES" ) : ( print OUT_info "\tNO" ) ;
				( $href_hoh_info_of_infile_of_type->{supermatrix}{convertedphy} == 1 ) ? ( print OUT_info "\tYES" ) : ( print OUT_info "\tNO" ) ;
				
				
				
				# (2)
				print OUT_info "\n\nFILE\tSTART POSITION SUPERMATRIX\tEND POSITION SUPERMATRIX\tN SEQUENCES\tSEQUENCE TYPE\tN CHARACTERS" ;
				
				for my $state ( "N MISSING DATA [?]", "N INDELS", "N AMBIGUITIES", "N PARSIMONY INFORMATIVE SITES", "N GC CONTENT", "N CODING CHARACTERS"	){ print OUT_info "\t", $state }
				for my $state ( qw/A C G T U N Y R W S K M D V H B F L I P Q E/												){ print OUT_info "\t", $state }
				
				for my $infile ( sort keys %$href_hol_seq_of_tax_of_infile ){
					
					( my $file = $infile ) =~ s/.fas$|.FASTA$|.aln$|.phy// ;
					
					print OUT_info	"\n", $file ;
					
					for my $state ( qw/seqstart seqend Ntaxa seqtype Nstates/			){ print OUT_info "\t", $href_hoh_info_of_infile_of_type->{$infile}{$state} }
					
					for my $state ( qw/Nmissingstate Ngapstate Nambstate infstate/		){
						
						( $href_hoh_info_of_infile_of_type->{$infile}{$state} ) ? ( print OUT_info "\t", $href_hoh_info_of_infile_of_type->{$infile}{$state} ) : ( print OUT_info "\t0" )
					}
					
					if ( $href_hoh_info_of_infile_of_type->{$infile}{seqtype} eq 'nu' ){
							
							( $href_hoh_info_of_infile_of_type->{$infile}{GC} ) ? ( print OUT_info "\t", $href_hoh_info_of_infile_of_type->{$infile}{GC} ) : ( print OUT_info "\t0" )
						}
					else{  print OUT_info "\tonly nuc data" }
					
					( $href_hoh_info_of_infile_of_type->{$infile}{Ninf_sites} ) ? ( print OUT_info "\t", $href_hoh_info_of_infile_of_type->{$infile}{Ninf_sites} ) : ( print OUT_info "\t0" ) ;
					
					
					for my $state ( qw/A C G T U N Y R W S K M D V H B F L I P Q E/		){
						
						( $href_hoh_info_of_infile_of_type->{$infile}{$state} ) ? ( print OUT_info "\t", $href_hoh_info_of_infile_of_type->{$infile}{$state} ) : ( print OUT_info "\t0" )
					}
				}
			}
			##############################
			
			
			
			##############################
			## If supermatrix is not printed out
			## (1) ...print infile name and sequence infos
			else{
				
				# (1)
				print OUT_info "File\tN SEQUENCES\tN SITES\tSEQUENCE TYPE\tN CHARACTERS" ;
				
				for my $state ( "N CODING CHARACTERS", "N MMISSING DATA [?]", "N INDELS", "N AMBIGUITIES", "N PARSIMONY INFORMATIVE SITES", "P GC CONTENT [%]"	){ print OUT_info "\t", $state }
				
				for my $state ( qw/A C G T U N Y R W S K M D V H B F L I P Q E/																				){ print OUT_info "\t", $state }
				
				for my $infile ( sort keys %$href_hol_seq_of_tax_of_infile ){
					
					( my $file = $infile ) =~ s/.fas$|.FASTA$|.aln$|.phy// ;
					
					print OUT_info	"\n", $file ;
					
					for my $state ( qw/Ntaxa seqlength seqtype Nstates/							){ print OUT_info "\t", $href_hoh_info_of_infile_of_type->{$infile}{$state} }
					
					for my $state ( qw/infstate Nmissingstate Ngapstate Nambstate/		){
						
						( $href_hoh_info_of_infile_of_type->{$infile}{$state} ) ? ( print OUT_info "\t", $href_hoh_info_of_infile_of_type->{$infile}{$state} ) : ( print OUT_info "\t0" )
					}
					
					( $href_hoh_info_of_infile_of_type->{$infile}{Ninf_sites} ) ? ( print OUT_info "\t", $href_hoh_info_of_infile_of_type->{$infile}{Ninf_sites} ) : ( print OUT_info "\t0" ) ;
					
					for my $state ( qw/GC/		){
						
						if ( $href_hoh_info_of_infile_of_type->{$infile}{seqtype} eq 'nu' ){
							
							if ( $href_hoh_info_of_infile_of_type->{$infile}{$state} ){
								
								print OUT_info "\t", ( sprintf "%.3f", ( $href_hoh_info_of_infile_of_type->{$infile}{$state} / $href_hoh_info_of_infile_of_type->{$infile}{Nstates} ) * 100 )
							}
							
							else{ print OUT_info "\t", 0 }
						}
						else{  print OUT_info "\tonly nuc data" }
					}
					
					for my $state ( qw/A C G T U N Y R W S K M D V H B F L I P Q E/	){
						
						( $href_hoh_info_of_infile_of_type->{$infile}{$state} ) ? (  print OUT_info "\t", $href_hoh_info_of_infile_of_type->{$infile}{$state} ) : ( print OUT_info "\t0" )
					}
				}
			}
			############################################################ END SEQUENCE INFO PRINT OUT
			
			
			
			##############################
			## If single data sets have been printed out in converted file format
			## print filename & output formats
			## Only for input files, not for supermatrix output
			if ( $aref_parameter_fil->[0] =~ /Convert/ ){
				
				print OUT_info "\n\nFile\tL SEQUENCE [bp]\tCONVERTED TO FASTA\tCONVERTED TO NEXUS\tCONVERTED TO PHYLIP" ;
				
				for my $infile ( sort keys %$href_hol_seq_of_tax_of_infile ){
					
					unless ( $infile eq 'supermatrix'){
						
						( my $file = $infile ) =~ s/.fas$|.aln$|.FASTA$|.phy$// ;
						
						print OUT_info "\n", $file, "\t", $href_hoh_info_of_infile_of_type->{$infile}{seqlength} ;
						
						( $href_hoh_info_of_infile_of_type->{$infile}{convertedfas} == 1 ) ? ( print OUT_info "\tYES" ) : ( print OUT_info "\tNO" ) ;
						( $href_hoh_info_of_infile_of_type->{$infile}{convertednex} == 1 ) ? ( print OUT_info "\tYES" ) : ( print OUT_info "\tNO" ) ;
						( $href_hoh_info_of_infile_of_type->{$infile}{convertedphy} == 1 ) ? ( print OUT_info "\tYES" ) : ( print OUT_info "\tNO" ) ;
					}
				}
			}
			##############################
			
			
			
			############################################################ START STRUCTURE INFO PRINT OUT
			## Print structure info if structure sequence found
			## (1) ...basal Structure info in common info file
			## (2) ...detailed structure info in separate structure info file
			## defined files for extra structure info output file are stored in @outfiles
			my @outfiles; 
			if ( $$sref_structure_seq ){
				
				
				############################################################ START SUPERMATRIX FILE HANDLING 
				## If File Handling option includes supermatrix output
				## (1) ...define supermatrix for increased structure info output in separate info file
				## (2) ...print basal structure info of supermatrix in info file
				if ( $aref_parameter_fil->[0] =~ /Supermatrix/ ){
					
					
					##############################
					## (1)
					## Define supermatrix print of extra structure info in separate outfile
					@outfiles = 'supermatrix' ;
					##############################
					
					
					
					##############################
					## (2) ...basal Structure info in common info file
					print OUT_info	"\n",
									"\n",
									"FASconCAT Supermatrix Structure Info\n",
									"Name Structure Sequence\tN Unpaired Positions (Loops)\tUnpaired Positions (Loops) [%]\tNumber Paired Positions (Stems)\tPaired Positions (Stems) [%]\n"
					;
					
					print OUT_info	$href_hoh_info_of_infile_of_type->{supermatrix}{taxstruct	}, "\t",
									$href_hoh_info_of_infile_of_type->{supermatrix}{Nloops		}, "\t",
									$href_hoh_info_of_infile_of_type->{supermatrix}{Ploops		}, "\t",
									$href_hoh_info_of_infile_of_type->{supermatrix}{Nstems		}, "\t",
									$href_hoh_info_of_infile_of_type->{supermatrix}{Pstems		}, "\n"
					;
					##############################
				}
				############################################################ ENND SUPERMATRIX FILE HANDLING 
				
				
				
				############################################################ START SINGLE INFILE HANDLING 
				## If File Handling option includes single infile output
				## (1) ...define all infiles for print of extra structure info in separate outfile
				## (2) ...print Structure info in common info file
				if ( $aref_parameter_fil->[0] =~ /Convert/ ){
					
					
					##############################
					## (1) ...define all infiles for print of extra structure info in separate outfile
					for my $infile ( sort {$a<=>$b} keys %$href_hol_seq_of_tax_of_infile ){ if ( $href_hoh_info_of_infile_of_type->{$infile}{taxstruct} ){ push @outfiles, $infile unless $infile eq 'supermatrix'} }
					##############################
					
					
					
					##############################
					## (2) ...print Structure info in common info file
					print OUT_info	"\n",
									"\n",
									"FASconCAT Supermatrix Structure Info Of Infiles Printed To ", $$sref_output_name_structure ,"\n"
					##############################
				}
				############################################################ END SINGLE INFILE HANDLING 
				
				
				
				############################################################ START OPEN OUT EXTRA STRUCTURE INFO FILE AND PRINT OUT
				## print detailed structure info in separate structure info file
				open  OUT_struc, ">$$sref_output_name_structure" ||  warn "\n\t!FILE-ERROR!: Cannot OPEN OUT ", $$sref_output_name_structure, "!\n" ;
				
				for my $infile ( sort {$a<=>$b} @outfiles ){
					
					
					##############################
					## Store identified sequences and associated taxon names of given infile in...
					## ...%sequence_of_taxon -> key: original taxon name; value: associated sequence
					my $aref_list_taxon_sequence_of_file	= exists( $href_hol_seq_of_tax_of_infile->{$infile}) ? \@{$href_hol_seq_of_tax_of_infile->{$infile}} :( ) ;
					my %sequence_of_taxon						= @$aref_list_taxon_sequence_of_file ;
					##############################
					
					
					
					##############################
					## Substitute file suffix
					( my $file = $infile ) =~ s/.fas$|.FASTA$|.aln$|.phy// ;
					##############################
					
					
					
					##############################
					## Print OUT extra structure info in separate outfile
					print OUT_struc	"FASconCAT Structure Info File:\t"	, $file, "\n",
										"Structure String\t"				, $sequence_of_taxon{$href_hoh_info_of_infile_of_type->{$infile}{taxstruct}}, "\n",
										"Loop positions\t"					, $href_hoh_info_of_infile_of_type->{$infile}{Loops		}, "\n",
										"Stem positions (paired)\t"			, $href_hoh_info_of_infile_of_type->{$infile}{Stems		}, "\n",
										"Stem positions (sorted)\t"			, $href_hoh_info_of_infile_of_type->{$infile}{Pairlist		}, "\n\n",
										"N Unpaired Positions (Loops):\t"	, $href_hoh_info_of_infile_of_type->{$infile}{Nloops		}, "\n",
										"Unpaired Positions (Loops) [%]\t"	, $href_hoh_info_of_infile_of_type->{$infile}{Ploops		}, "\n",
										"Number Paired Positions (Stems)\t"	, $href_hoh_info_of_infile_of_type->{$infile}{Nstems		}, "\n",
										"Paired Positions (Stems) [%]\t"	, $href_hoh_info_of_infile_of_type->{$infile}{Pstems		}, "\n\n"
					;
					##############################
				}
				
				##############################
				## Close extra structure info outfile
				close OUT_struc ;
				##############################
				
				############################################################ END OPEN OUT EXTRA STRUCTURE INFO FILE AND PRINT OUT
			}
			############################################################ END SUPERMATRIX FILE HANDLING 
			
			
			
			##############################
			## Print... 
			## (1) ...N of taxon concatenations (total)
			## (2) ...N informative taxon sequences
			## (3) ...N missing taxon sequences
			if ( $aref_parameter_fil->[0] =~ /Supermatrix/ ){
				print OUT_info "\n",
								"\nTaxon check\n",
								"Taxon\tN CONCATENATED SEQs\tN INFORMATIVE SEQs\tN MISSING SEQs\n"
				;
				
				my %supermatrix_seq_of_taxon = @{$href_hol_seq_of_tax_of_infile->{supermatrix}} ;
				for my $taxon ( sort {$a<=>$b} keys %supermatrix_seq_of_taxon ){
					
					unless ( $href_hoh_info_of_infile_of_type->{supermatrix}{Nseqconcat	}{$taxon} ){ $href_hoh_info_of_infile_of_type->{supermatrix}{Nseqconcat	}{$taxon} = 0 }
					unless ( $href_hoh_info_of_infile_of_type->{supermatrix}{Nsequnconcat	}{$taxon} ){ $href_hoh_info_of_infile_of_type->{supermatrix}{Nsequnconcat	}{$taxon} = 0 }
					
					
					print OUT_info $taxon,
												"\t", $href_hoh_info_of_infile_of_type->{supermatrix}{taxconcat		}{$taxon},	# (1)
												"\t", $href_hoh_info_of_infile_of_type->{supermatrix}{Nseqconcat	}{$taxon},	# (2)
												"\t", $href_hoh_info_of_infile_of_type->{supermatrix}{Nsequnconcat	}{$taxon},	# (3)
												"\n"
				}
			}
			##############################
			
			
			
			##############################
			## Terminal Info Print...
			## (1) ...Site Ranges of given infiles in the concatenated supermatrix (if any...)
			## (2) ...% Missing site states (X or ?) of given infiles and... 
			## (3) ...the concatenated supermatrix (if given)
			
			# (1)
			if ( $aref_parameter_fil->[0] =~ /Supermatrix/ ){
				
				print "\n\n\tSupermatrix Site Range Of File...\n\t" ;
				for my $infile ( sort keys %$href_hol_seq_of_tax_of_infile ){
					
					( my $file = $infile ) =~ s/.fas$|.aln$|.FASTA$|.phy$// ;
					
					unless ( $file eq 'supermatrix' ){ print $file, ":\tSite Pos.", $href_hoh_info_of_infile_of_type->{$infile}{seqstart}, " to ", $href_hoh_info_of_infile_of_type->{$infile}{seqend}, "\n\t" }
				}
			}
			
			# (2)
			print "\n\tPercent Missing Data (X or ?) Of File..." ;
			for my $infile ( sort keys %$href_hol_seq_of_tax_of_infile ){
				
				unless ( $infile eq 'supermatrix' ){
					
					( my $file = $infile ) =~ s/.fas$|.FASTA$|.aln$|.phy// ;
					
					if ( $href_hoh_info_of_infile_of_type->{$infile}{Nmissingstate} ){
						
						print "\n\t", $file, ":\t", ( sprintf "%.3f", ( $href_hoh_info_of_infile_of_type->{$infile}{Nmissingstate} / $href_hoh_info_of_infile_of_type->{$infile}{Nstates} ) * 100 ), "%"
					}
					
					else{ print "\n\t", $file, ":\t", 0, "%" }
				}
			}
			
			# (3)
			if ( $aref_parameter_fil->[0] =~ /Supermatrix/ ){
				
				if ( $href_hoh_info_of_infile_of_type->{supermatrix}{Nmissingstate} ){
					
					print "\n\tSupermatrix:\t", ( sprintf "%.3f", ( $href_hoh_info_of_infile_of_type->{supermatrix}{Nmissingstate} / $href_hoh_info_of_infile_of_type->{supermatrix}{Nstates} ) * 100 ), "%"
				}
				
				else{ print "\n\tSupermatrix:\t", 0, "%" }
			}
			##############################
			
			
			
			##############################
			# Print info header
			print OUT_info "\n"							,
							"\n------------------------------------------------------------" ,
							"\nDefined Options:"			,
							"\n-----------------"			,
							"\nREAD IN ALL files:\t"		,	$aref_parameter_all->[0],
							"\nREAD IN SINGLE files:\t"		,	$aref_parameter_all->[1],
							"\nFILE HANDLING:\t"			,	$aref_parameter_fil->[0],
							"\nSEQUENCE TRANSLATION\t"		,	$aref_parameter_tra->[0],
							"\nEXCLUDING 3rd POSITION\t"	,	$aref_parameter_3rd->[0],
							"\nRY CODING\t"					,	$aref_parameter_ryc->[0],
							"\nBUILD CONSENSUS\t"			,	$aref_parameter_con->[0],
							"\nRENAME SEQUENCE NAMES\t"		,	$aref_parameter_ren->[0],
							"\nABSENT SEQUENCE CODING\t"	,	$aref_parameter_mis->[0],
							"\n"							,
							"\n-----------------"			,
							"\nOUTPUT FORMAT:"				,
							"\nPHYLIP:\t"					,	$aref_parameter_phy->[0],
							"\nNEXUS:\t"					,	$aref_parameter_nex->[0],
							"\nFASTA:\t"					,	$aref_parameter_fas->[0],
							"\nPARTITION:\t"				,	$aref_parameter_prt->[0],
							"\nPARSIMONY:\t"				,	$aref_parameter_par->[0],
							"\nPROTTEST:\t"					,	$aref_parameter_pro->[0],
							"\n------------------------------------------------------------\n"
			;
			##############################
			
			
			
			##############################
			## close info file
			close OUT_info ;
			##############################
		}
		
		sub fas_out{
			
			my $aref_outfile_names				= $_[0] ;	# list of filenames for print OUT														-> IN (defined) / OUT (unchanged)
			my $href_hol_seq_of_tax_of_infile	= $_[1] ;	# key: input file name; value ; list of taxa and corresponding sequences				-> IN (defined) / OUT (unchanged)
			my $href_hoh_info_of_infile_of_type	= $_[2] ;	# key1: infilename, key2: type of info (e.g. 'seqlength'); info type associated value	-> IN (defined) / OUT (changed)
			my $sref_parsimony_code				= $_[3] ;	# if 'parsimony' do not set conversion counter +1										-> IN (defined) / OUT (unchanged)
			
			
			############################## START FASTA PRINT OUT OF DEFINDED FILE DATA
			## for each defined outfile
			## extract taxon and sequence info from %$href_hol_seq_of_tax_of_infile
			## open out fasta file of infile and print out fasta formated taxon and sequence data
			for my $file ( @$aref_outfile_names ){
				
				unless ( $$sref_parsimony_code eq 'parsimony' ){ $href_hoh_info_of_infile_of_type->{$file}{convertedfas} = 1 }
				
				my @data = exists($href_hol_seq_of_tax_of_infile->{$file}) ? @{$href_hol_seq_of_tax_of_infile->{$file}} :( ) ;
				
				( my	$outfile = $file ) =~ s/.phy$|.aln$|.FASTA$|.fas$// ;
						$outfile = 'FcC_'.$outfile.".fas" ;
				
				open  OUT_fas,	">$outfile" || warn "\n\t!FILE-ERROR!: Cannot OPEN OUT ", $outfile, "!\n" ;
				while ( @data ){
					
					my $taxon		= shift @data ;
					my $sequence	= shift @data ;
					
					print OUT_fas ">", $taxon, "\n", $sequence, "\n"
				}
				close OUT_smatrix
			}
			############################## END FASTA PRINT OUT OF DEFINDED FILE DATA
		}
		
		sub nexus_out{
			
			my $aref_outfile_names					= $_[0] ;	# list of filenames for print OUT														-> IN (defined) / OUT (unchanged)
			my $sref_nexus_setting_out				= $_[1] ;	# nexus output parameter ('BLOCK' or 'MrBAYES')											-> IN (defined) / OUT (unchanged)
			my $href_hol_seq_of_tax_of_infile		= $_[2] ;	# key: input file name; value ; list of taxa and corresponding sequences				-> IN (defined) / OUT (unchanged)
			my $href_hoh_info_of_infile_of_type		= $_[3] ;	# key1: infilename, key2: type of info (e.g. 'seqlength'); info type associated value	-> IN (defined) / OUT (unchanged)
			my $sref_parsimony_code					= $_[4] ;	# if 'parsimony' do not set conversion counter +1										-> IN (defined) / OUT (unchanged)
			
			
			# MrBAYES - PARAMETER
			my $partition_looms	= "partition looms = 2: loops, stems;" ;
			my $set_partition		= "set partition = looms;" ;
			my $lset_applyto1		= "lset applyto= (1) nucmodel= 4by4;" ;
			my $lset_applyto2		= "lset applyto= (2) nucmodel= dublet;" ;
			my $lset_applyto12	= "lset applyto= (1,2) nst= 6 rates= gamma;" ;
			
			
			
			
			# NEXUS - BLOCK - PARAMETER
			my $format			= "format datatype=dna interleave" ;
			my $missing		= "missing=-;" ;
			my $autoclose		= "set autoclose= yes" ;
			my $nowarn			= "nowarn=yes;" ;
			my $state_freq	= "unlink statefreq= (all)" ;
			my $tratio			= "tratio= (all)" ;
			my $shape			= "shape= (all);" ;
			my $mcmc_gen		= "mcmc ngen= 2000" ;
			my $print_freq	= "printfreq= 100" ;
			my $sample_freq	= "samplefreq= 100" ;
			my $nchains		= "nchains= 4" ;
			my $save_brlens	= "savebrlens= yes" ;
			my $sump_burnin	= "sump burnin= 20" ;
			my $sumt_burnin	= "sumt burnin= 10" ;
			my $nruns			= "nruns= 2;" ;
			
			
			
			
			############################################################ START NEXUS PRINT OUT OF DEFINDED FILE DATA
			## for each defined outfile
			## extract taxon and sequence info from %$href_hol_seq_of_tax_of_infile
			## open out phylip file of infile and print out phylip formated taxon 
			## and sequence data according to chosen phylip format
			for my $file ( @$aref_outfile_names ){
				
				
				unless ( $$sref_parsimony_code eq 'parsimony' ){ $href_hoh_info_of_infile_of_type->{$file}{convertednex} = 1 }
				
				my @data = exists($href_hol_seq_of_tax_of_infile->{$file}) ? @{$href_hol_seq_of_tax_of_infile->{$file}} :( ) ;
				
				
				
				##############################
				## (1) ...assign list of taxa and corresponding sequences to hash %sequence_of_taxon
				## (2) ...determine number of sequence positions 
				## (3) ...determine number of sequences
				## (4) ...if nexus is given out as 'MrBAYES' outputfile and structure sequence found, decrease number of sequences by -1
				my %sequence_of_taxon	= @data										;	# (1)
				my $sequence_length	= length $sequence_of_taxon{$data[0]}	;	# (2)
				my $Ntaxa				= keys %sequence_of_taxon					;	# (3)
				
				# (4)
				if ( ( $href_hoh_info_of_infile_of_type->{$file}{taxstruct} ) && ( $$sref_nexus_setting_out eq 'MrBAYES' ) ) { $Ntaxa-- }
				##############################
				
				
				
				############################################################ START GENERATING NEXUS SEQUENCE BLOCKS (INTERLEAVED)
				## Generate Nexus sequence blocks interleaved formated
				## always 20 base pairs / block & 5 blocks in a row
				##
				## -----------------------
				## e.g.:
				## 			BLOCK_1 (1-20) | BLOCK_2 (21-40) | BLOCK_3 (41-60) | BLOCK_4 (61-80) | BLOCK_5 (81-100)
				## tax_1
				## tax_2
				## tax_3
				## tax_n
				## 			BLOCK_1 (101-120) | BLOCK_1 (121-140) | BLOCK_3 (141-160) | BLOCK_4 (161-180) | BLOCK_5 (181-200)
				## tax_1
				## tax_2
				## tax_3
				## tax_n
				## -----------------------
				##
				## For each taxon store...
				## (1) ...pack blocks of 20 sequence positions in @block_rows_nex
				## (2) ...join always 5 blocks of @block_rows_nex to one string in %hoh_seq_parts_nex -> key1: taxon; key2: counter of joined blocks; value: joined blocks
				my %hoh_seq_parts_nex ; my $counter_blocks ;
				for my $nex_tax ( sort {$a<=>$b} keys %sequence_of_taxon ){
					
					
					my (
							@block_rows_nex,	# stores sequence blocks of 20 bp
					) ;
					
					$counter_blocks = 0 ;
					
					# (1)
					for ( my $j=0; $j<= $sequence_length; $j+=20 ){ push @block_rows_nex, substr($sequence_of_taxon{$nex_tax}, $j, 20) }
					
					# (2)
					for ( my $k=1; $k<= @block_rows_nex; $k+=5  ){
						
						$counter_blocks++ ;
						$hoh_seq_parts_nex{$nex_tax}{$counter_blocks} = ( join " ", ( $block_rows_nex[$k-1], $block_rows_nex[$k], $block_rows_nex[$k+1], $block_rows_nex[$k+2], $block_rows_nex[$k+3] ) ) ;
					}
				}
				############################################################ END GENERATING NEXUS SEQUENCE BLOCKS (INTERLEAVED)
				
				
				
				############################################################ START OPEN OUTFILE & PRINT OUT NEXUS HEADER
				## Open OUT nexus output file
				## Print NEXUS header
				( my	$outfile = $file ) =~ s/.fas$|.aln$|.FASTA$|.phy$// ;
						$outfile = 'FcC_'.$outfile.".nex" ;
					
				open  OUT_nex,	">$outfile" || warn "\n\t!FILE-ERROR!: Cannot OPEN OUT ", $outfile, "!\n" ;
				print OUT_nex	
								"#NEXUS\n",
								"\n",
								"begin data;\n",
								"dimensions ntax= ", $Ntaxa, " nchar= ", $sequence_length, ";\n" ,
								$format, " ", $missing, "\n",
								"matrix\n"
				;
				############################################################ END OPEN OUTFILE & PRINT OUT NEXUS HEADER
				
				
				
				############################################################ START PRINT OUT NEXUS BLOCKS
				## Print NEXUS sequence blocks
				## per row -> 5 blocks with 20 characters each
				for my $row ( 1 .. $counter_blocks ){
					
					for my $taxon ( sort {$a<=>$b} keys %hoh_seq_parts_nex ){
						
						unless ( $taxon eq $href_hoh_info_of_infile_of_type->{$file}{taxstruct} ){ print OUT_nex $taxon, "  ",  $hoh_seq_parts_nex{$taxon}{$row}, "\n" }
					}
					
					print OUT_nex "\n" ;
				}
				############################################################ END PRINT OUT NEXUS BLOCKS
				
				
				
				############################################################ START PRINT OUT NEXUS COMMAND BELOW BLOCK
				## Print NEXUS commands below nexus block for 'BLOCK' output
				( my $log_outfile = $outfile ) =~s/.nex$/.txt/ ;
				if ( $$sref_nexus_setting_out eq 'BLOCK' ){
					
					print OUT_nex 
									";\n",
									"end;\n",
									"log start file= ", $log_outfile, " append;\n",
									"end ;\n"
				}
				##############################
				
				
				
				##############################
				## Print NEXUS commands below nexus block for 'MrBAYES' output
				elsif ( $$sref_nexus_setting_out eq 'MrBAYES' ){
					
					my $bayes_block_string ;
					
					if   ( $href_hoh_info_of_infile_of_type->{$file}{taxstruct} ) {
						
							$bayes_block_string = "\n\nbegin mrbayes;\ncharset loops= ".$href_hoh_info_of_infile_of_type->{$file}{Loops}.";\ncharset stems= ".$href_hoh_info_of_infile_of_type->{$file}{Stems}.";\n\npairs ".$href_hoh_info_of_infile_of_type->{$file}{Pairs}.";\n\n".$partition_looms."\n\n".$set_partition."\n".$lset_applyto1."\n".$lset_applyto2."\n".$lset_applyto12 ;
					}
					
					else {	$bayes_block_string = "\n\nbegin mrbayes;\n lset nst= 6 rates= gamma;\n\n" }
					
					print OUT_nex 
									";\n",
									"end;\n",
									"log start file= ", $log_outfile, " append;\n",
									"\n",
									$bayes_block_string, "\n",
									"showmodel;\n",
									"$autoclose ", $nowarn, "\n",
									$state_freq, $tratio, " ", $shape, "\n",
									$mcmc_gen, " ", $print_freq, " ", $sample_freq, " ", $nchains, " ", $save_brlens, " filename= FcC_mrbayes_output.txt ;\n",
									$sump_burnin, " ", $nruns, "\n",
									$sumt_burnin, " ", $nruns, "\n",
									"\n",
									"quit;\n"
					;
				}
				##############################
				
				
				
				##############################
				## BUG REPORT
				else { warn "\n\tBUG-ERROR: Cannot assign NEXUS output parameter ", $$sref_nexus_setting_out, "!\n"  }
				##############################
				
				############################################################ END PRINT OUT NEXUS COMMAND BELOW BLOCK
				
				
				
				##############################
				## Close NEXUS outfile
				close OUT_nex
				##############################
			}
			############################################################ END NEXUS PRINT OUT OF DEFINDED FILE DATA
		}
		
		sub nexus_out_prt{
			
			my $sref_nexus_setting_out				= $_[0] ;	# nexus output parameter ('BLOCK' or 'MrBAYES')											-> IN (defined) / OUT (unchanged)
			my $sref_partition_setting				= $_[1] ;	# nexus output parameter ('Supermatrix' or '3rd')										-> IN (defined) / OUT (unchanged)
			my $href_hol_seq_of_tax_of_infile		= $_[2] ;	# key: input file name; value ; list of taxa and corresponding sequences				-> IN (defined) / OUT (unchanged)
			my $href_hoh_info_of_infile_of_type	= $_[3] ;	# key1: infilename, key2: type of info (e.g. 'seqlength'); info type associated value	-> IN (defined) / OUT (unchanged)
			
			if ( $$sref_partition_setting eq 'Supermatrix' ){
				
				##############################
				# OPEN OUT nexus partition file
				open OUTnex, ">FcC_supermatrix_partition.nex" || warn "\n\t!FILE-ERROR!: Cannot OPEN OUT  FcC_supermatrix_partition.nex!\n" ;
				##############################
				
				
				
				##############################
				# PRINT HEADER
				print OUTnex	"#NEXUS\n",
								"BEGIN DATA;\n",
								"\tDIMENSIONS NTAX=", $href_hoh_info_of_infile_of_type->{supermatrix}{Ntaxa}, " NCHAR=", $href_hoh_info_of_infile_of_type->{supermatrix}{seqlength}, ";\n"
				;
				
				if ( $href_hoh_info_of_infile_of_type->{supermatrix}{seqtype} eq 'nu' )	{ print OUTnex	"\tFORMAT DATATYPE = DNA GAP = - MISSING = ? INTERLEAVE;\n" }
				else 																			{ print OUTnex	"\tFORMAT DATATYPE = Protein GAP = - MISSING = ? INTERLEAVE;\n" }
				
				print OUTnex	"\tMATRIX\n";
				##############################
				
				
				
				##############################
				## Store identified sequences and associated taxon names of given infile in...
				## ...%sequence_of_taxon -> key: original taxon name; value: associated sequence
				my $aref_list_taxon_sequence_of_file	= exists( $href_hol_seq_of_tax_of_infile->{supermatrix}) ? \@{$href_hol_seq_of_tax_of_infile->{supermatrix}} :( ) ;
				my %sequence_of_taxon						= @$aref_list_taxon_sequence_of_file ;
				##############################
				
				
				
				##############################
				# PRINT PARTITION SEQUENCE BLOCKS
				for my $infile ( sort {$a<=>$b} keys %$href_hol_seq_of_tax_of_infile ){
					
					unless ( $infile eq 'supermatrix' ){
						
						my $startposition	= $href_hoh_info_of_infile_of_type->{$infile}{seqstart} ;
						my $endposition	= $href_hoh_info_of_infile_of_type->{$infile}{seqend} ;
						my $length			= $href_hoh_info_of_infile_of_type->{$infile}{seqlength} ;
						#my $Nchar			= $href_hoh_info_of_infile_of_type->{supermatrix}{Ntaxa} * $length ;
						
						print OUTnex "[",$infile, "\tdimensions nchar=", $length, "]\n" ;
						
						for my $taxon ( sort keys %sequence_of_taxon ){
							
							my $seq_substring = substr( $sequence_of_taxon{$taxon}, $startposition-1, $length ) ;
							
							print OUTnex	$taxon, "\t\t", $seq_substring, "\n"
						}
						print OUTnex	"\n" ;
					}
				}
				print OUTnex	";\n\nEND;\n" ;
				##############################
				
				
				
				##############################
				# PRINT MrBayes COMMAND BLOCKS
				if ( $$sref_nexus_setting_out eq 'MrBAYES' ){
					
					##############################
					# print MrBayes header
					print OUTnex	"\n\nBEGIN mrbayes;\n",
									"\tset autoclose=yes nowarn=yes;\n",
									"\tlog start filename=FcC_supermatrix_partition_log.txt;\t\t\t\t\t[log output]\n",
									"\t[outgroup outgroupname;]\t\t\t\t\t\t\t\t\t\t\t\t\t[define outgroup]\n"
					;
					##############################
					
					
					
					##############################
					# print defined char sets
					my @charsets;
					for my $infile ( sort {$a<=>$b} keys %$href_hol_seq_of_tax_of_infile ){
						
						unless ( $infile eq 'supermatrix' ){
							
							print OUTnex "\tcharset ", $infile, " = ", $href_hoh_info_of_infile_of_type->{$infile}{seqstart}, "-", $href_hoh_info_of_infile_of_type->{$infile}{seqend}, ";\n" ;
							push @charsets, $infile ;
						}
					}
					##############################
					
					
					
					##############################
					# print mrbayes setup
					my $Ncharsets	= @charsets ;
					my $charsets	= join ", ", @charsets ;
					
					print OUTnex	"\tpartition genes = ", $Ncharsets, ": ", $charsets, ";\t\t[define a partition scheme]\n",
									"\tset partition=genes;\t\t\t\t\t\t\t\t\t\t\t[set the partition scheme]\n",
									"\tlset applyto=(all) nst=6 rates=invgamma;\t\t\t\t\t\t[set a GTR+I+G model for all partitions]\n",
									"\tunlink shape=(all) pinvar=(all) statefreq=(all) revmat=(all);\t[specify mixed model--how parameters are shared across partitions...]\n",
									"\tprset applyto=(all) ratepr=variable;\t\t\t\t\t\t\t[...and allow rate to vary across partitions]\n",
									"\tmcmc ngen=10000000 printfreq=1000 samplefreq=1000\t\t\t\t[specify MCMC details...run length and thinning]\n",
									"\tdiagnfreq=1000 diagnstat=maxstddev\t\t\t\t\t\t\t\t[...MCMC diagnostics...]\n",
									"\tcheckpoint=yes checkfreq=100000\t\t\t\t\t\t\t\t\t[...checkpointing so runs can be restarted...]\n",
									"\tnchains=4 temp=0.1\t\t\t\t\t\t\t\t\t\t\t\t[specify MCMCMC details]\n",
									"\tsavebrlens=yes\n",
									"\tfilename=FcC_supermatrix_partition;\n",
									"\tsumt filename=FcC_supermatrix_partition;\n",
									"\tsump filename=FcC_supermatrix_partition;\n",
									"END;\n\n"
					;
					##############################
				}
				close OUTnex;
				##############################
			}
		}
		
		sub phylip_out{
			
			my $aref_outfile_names				= $_[0] ;	# list of filenames for print OUT														-> IN (defined) / OUT (unchanged)
			my $href_hol_seq_of_tax_of_infile	= $_[1] ;	# key: input file name; value ; list of taxa and corresponding sequences				-> IN (defined) / OUT (unchanged)
			my $sref_phylip_setting_out			= $_[2] ;	# phylip output parameter ('strict' or 'relaxed')										-> IN (defined) / OUT (unchanged)
			my $href_hoh_info_of_infile_of_type	= $_[3] ;	# key1: infilename, key2: type of info (e.g. 'seqlength'); info type associated value	-> IN (defined) / OUT (changed)
			my $sref_parsimony_code				= $_[4] ;	# if 'parsimony' do not set conversion counter +1										-> IN (defined) / OUT (unchanged)
			
			
			############################################################ START PHYLIP PRINT OUT OF DEFINDED FILE DATA
			## for each defined outfile
			## extract taxon and sequence info from %$href_hol_seq_of_tax_of_infile
			## open out phylip file of infile and print out phylip formated taxon 
			## and sequence data according to chosen phylip format
			for my $file ( @$aref_outfile_names ){
				
				unless ( $$sref_parsimony_code eq 'parsimony' ){ $href_hoh_info_of_infile_of_type->{$file}{convertedphy} = 1 }
				
				my @data = exists($href_hol_seq_of_tax_of_infile->{$file}) ? @{$href_hol_seq_of_tax_of_infile->{$file}} :( ) ;
				
				( my	$outfile = $file ) =~ s/.fas$|.aln$|.FASTA$// ;
						$outfile = 'FcC_'.$outfile.".phy" ;
				
				##############################
				## Open OUT phylip output file
				## Print phylip header
				## N taxa and N character states
				## ---------------
				## e.g.
				## 12 4356
				## ---------------
				my $Ntaxa	= @data / 2 ;
				my $Nchar	= length $data[1] ;
				
				open  OUTphy,	">$outfile" || warn "\n\t!FILE-ERROR!: Cannot OPEN OUT ", $outfile, "!\n" ;
				print OUTphy	"$Ntaxa $Nchar\n" ;
				##############################
				
				
				
				############################################################ START RELAXED PRINT OUT
				## PRINT OUT supermatrix RELAXED formated
				## taxon names are not restricted to 10 signs 
				## ---------------
				## e.g.
				## 12 4356
				## ---------------
				if ( $$sref_phylip_setting_out eq 'RELAXED' ){ 
					
					while ( @data ){
						
						my $taxon				= shift @data ;
						my $sequence_of_taxon	= shift @data ;
						
						( my $taxon_phy = $taxon ) =~ s/ /_/g ;
						
						print OUTphy $taxon_phy, "   ", $sequence_of_taxon, "\n"
					}
				}
				############################################################ END RELAXED PRINT OUT
				
				
				
				############################################################ START STRICT PHYLIP PRINT OUT
				## PRINT OUT supermatrix STRICT formated
				## taxon names are restricted to 10 signs
				## ---------------
				## e.g.
				## taxon_1   ACCCGGGTTTTAA...
				## taxon_2   ACCCGGTTTTTAA...
				## taxon_n   ACCCGGGTGGTAA...
				## ---------------
				elsif ( $$sref_phylip_setting_out eq 'STRICT' ){
					
					
					##############################
					## define counter of taxon names which have been restricted to 10 signs
					## If restricted taxon name equal previously restricted taxon name...
					## (1) ...set taxon counter of restricted taxon plus 1
					## (2) ...replace last signs of restricted taxon name in length of taxon counter number with counter number
					my %seen_taxon_restricted ;
					##############################
					
					
					
					############################################################ START TAXON & SEQUENCE PRINT OUT
					## while list of taxon names and sequences...
					## (1) ...assign taxon name
					## (2) ...determine number of signs in taxon name
					## (3) ...if number of taxon signs are greater -> 10...
					##		(3a) ...reduce taxon signs to 10
					while ( @data ){
						
						my $taxon					= shift @data	;	# (1)
						my $Nsigns_taxon			= length $taxon ;	# (2)
						
						############################################################ START HANDLING OF TAXON NAMES GREATR 10 SIGNS
						## (3) if number of taxon signs are greater -> 10...
						if ( $Nsigns_taxon > 10 ){
							
							
							
							##############################
							## (3a) reduce taxon signs to 10
							my $taxon_restricted = substr($taxon, 0, 10) ;
							$taxon_restricted =~ s/ /_/g ;
							##############################
							
							
							
							##############################
							## take taxon corresponding sequence from list @data
							my $sequence_of_taxon	= shift @data ;
							##############################
							
							
							
							##############################
							# PRINT OUT taxon and corresponding sequence
							print OUTphy $taxon_restricted, "   ", $sequence_of_taxon, "\n"
							##############################
						}
						############################################################ END HANDLING OF TAXON NAMES GREATR 10 SIGNS
						
						
						else{
							
							##############################
							## take taxon corresponding sequence from list @data
							my $sequence_of_taxon	= shift @data ;
							##############################
							
							
							
							##############################
							# PRINT OUT taxon and corresponding sequence
							print OUTphy $taxon, "   ", $sequence_of_taxon, "\n"
							##############################
						}
					}
					############################################################ END TAXON & SEQUENCE PRINT OUT
				}
				############################################################ END STRICT PHYLIP PRINT OUT
				
				
				
				##############################
				## CLOSE phylip output file
				close OUTphy;
				##############################
			}
			############################################################ END PHYLIP PRINT OUT OF DEFINDED FILE DATA
		}
		
		sub raxml_out_prt{
			
			my $sref_partition_setting				= $_[0] ;	# nexus output parameter ('Supermatrix' or '3rd')										-> IN (defined) / OUT (unchanged)
			my $href_hol_seq_of_tax_of_infile		= $_[1] ;	# key: input file name; value ; list of taxa and corresponding sequences				-> IN (defined) / OUT (unchanged)
			my $href_hoh_info_of_infile_of_type		= $_[2] ;	# key1: infilename, key2: type of info (e.g. 'seqlength'); info type associated value	-> IN (defined) / OUT (unchanged)
			my $sref_prottest_option				= $_[3] ;	# parameter for prottest analyses														-> IN (defined) / OUT (unchanged)
			
			
			############################################################ START SUPERMATRIX PARTITION FILE
			if ( $$sref_partition_setting eq 'Supermatrix' ){
				
				##############################
				# OPEN OUT fasta partition file
				open OUTfas, ">FcC_supermatrix_partition.txt" || warn "\n\t!FILE-ERROR!: Cannot OPEN OUT  FcC_supermatrix_partition.txt!\n" ;
				##############################
				
				
				
				##############################
				# print raxml partition file
				# if prottests have not been performed use default models for raxml partition file
				# otherwise use best bic models of prottests
				my $seqtype ;
				for my $infile ( sort keys %$href_hol_seq_of_tax_of_infile ){
					
					unless ( $infile eq 'supermatrix' ){
						
						( my $file = $infile ) =~ s/.phy$|.aln$|.fas$|.FASTA$// ;
						
						unless ( $$sref_prottest_option eq 'Supermatrix' ){ if ( $href_hoh_info_of_infile_of_type->{supermatrix}{seqtype} eq 'nu' ){ $seqtype = 'DNA' } else { $seqtype = 'LG' } }
						else{ if ( $href_hoh_info_of_infile_of_type->{$infile}{protmodel} ){ $seqtype = $href_hoh_info_of_infile_of_type->{$infile}{protmodel} } else { $seqtype = 'LG'; print "\n\tProtTest not executed correctly for ", $file, "!\n\tSet partition model to LG" } }
						
						print OUTfas $seqtype, ", ", $file, " = ", $href_hoh_info_of_infile_of_type->{$infile}{seqstart}, "-", $href_hoh_info_of_infile_of_type->{$infile}{seqend}, "\n"
					}
				}
				##############################
			}
			############################################################ END SUPERMATRIX PARTITION FILE
			
			
			
			##############################
			# close FASTA partition file
			close OUTfas
			##############################
		}
		
		sub parsimony_print{
			
			my $aref_outfile_names				= $_[0] ;	# list of filenames for print OUT														-> IN (defined) / OUT (unchanged)
			my $href_hol_seq_of_tax_of_infile	= $_[1] ;	# key: input file name; value ; list of taxa and corresponding sequences				-> IN (defined) / OUT (unchanged)
			my $href_hoh_info_of_infile_of_type	= $_[2] ;	# key1: infilename, key2: type of info (e.g. 'seqlength'); info type associated value	-> IN (defined) / OUT (unchanged)
			my $sref_parameter_fas				= $_[3] ;	# defined parameter option of FASTA print OUT											-> IN (defined) / OUT (unchanged)
			my $sref_parameter_nex				= $_[4] ;	# defined parameter option of NEXUS print OUT											-> IN (defined) / OUT (unchanged)
			my $sref_parameter_phy				= $_[5] ;	# defined parameter option of PHYLIP print OUT											-> IN (defined) / OUT (unchanged)
			
			
			##############################
			## new variabels for subroutine print out
			my %hol_inf_site_seq_of_tax_of_infile ;
			my @parsimony_file_names ;
			##############################
			
			
			
			############################################################ START EXTRACTION OF PARSIMONY INFORMATIVE SITE FOR EACH INFILE
			for my $file ( @$aref_outfile_names ){
				
				
				##############################
				## define parsimony outfile name for each file and store
				## outfilenames as list @parsimony_file_names
				my $parsimony_file = "parsim_inf_sites_".$file ;
				push @parsimony_file_names, $parsimony_file ;
				##############################
				
				
				
				##############################
				## count number of identified informative sites in %seen_inf_pos
				## of given infile
				my $aref_inf_sites = exists ($href_hoh_info_of_infile_of_type->{$file}{list_inf_pos}) ? \@{$href_hoh_info_of_infile_of_type->{$file}{list_inf_pos}} :( ) ; 
				my %seen_inf_pos ;
				for my $inf_pos ( @$aref_inf_sites ){ $seen_inf_pos{$inf_pos}++ }
				##############################
				
				
				
				##############################
				## generate for each taxon sequence
				## a new sequence string with only informative sites
				my @data = exists($href_hol_seq_of_tax_of_infile->{$file}) ? @{$href_hol_seq_of_tax_of_infile->{$file}} :( ) ;
				
				while ( @data ){
					
					my $taxon		= shift @data ;
					my $sequence	= shift @data ;
					
					my @sites = split "", $sequence ;
					$sequence = () ;
					for ( 0 .. @sites-1 ){ if ( $seen_inf_pos{$_} ){ $sequence .= $sites[$_] } }
					
					push @{$hol_inf_site_seq_of_tax_of_infile{$parsimony_file}}, ( $taxon, $sequence )
				}
				##############################
			}
			############################################################ END EXTRACTION OF PARSIMONY INFORMATIVE SITE FOR EACH INFILE
			
			
			
			############################################################ START PRINT OUT OF PARSIMONY INFORMATIVE SITE FOR EACH INFILE
			
			##############################
			## OUTPUT -> FASTA format 
			unless ( $$sref_parameter_fas eq 'NO' ){
				
				&fas_out(
							\@parsimony_file_names,					# list of filenames for print OUT														-> IN (defined) / OUT (unchanged)
							\%hol_inf_site_seq_of_tax_of_infile,	# key: input file name; value ; list of taxa and corresponding sequences				-> IN (defined) / OUT (unchanged)
							\%$href_hoh_info_of_infile_of_type,		# key1: infilename, key2: type of info (e.g. 'seqlength'); info type associated value	-> IN (defined) / OUT (unchanged)
							\'parsimony'							# if 'parsimony' do not set conversion counter +1										-> IN (defined) / OUT (unchanged)
				)
			}
			##############################
			
			
			
			##############################
			## Additional output (NEXUS 'BLOCK' or 'MrBAYES')
			unless ( $$sref_parameter_nex eq 'NO' ){
				
				&nexus_out(
							\@parsimony_file_names,					# list of filenames for print OUT														-> IN (defined) / OUT (unchanged)
							\$$sref_parameter_nex,					# nexus output parameter ('BLOCK' or 'MrBAYES')											-> IN (defined) / OUT (unchanged)
							\%hol_inf_site_seq_of_tax_of_infile,	# key: input file name; value ; list of taxa and corresponding sequences				-> IN (defined) / OUT (unchanged)
							\%$href_hoh_info_of_infile_of_type,		# key1: infilename, key2: type of info (e.g. 'seqlength'); info type associated value	-> IN (defined) / OUT (unchanged)
							\'parsimony'							# if 'parsimony' do not set conversion counter +1										-> IN (defined) / OUT (unchanged)
				)
			}
			##############################
			
			
			##############################
			## Additional output (phylip STRICT or RELAXED)
			unless ( $$sref_parameter_phy eq 'NO' ){
				
				&phylip_out(
							\@parsimony_file_names,					# list of filenames for print OUT														-> IN (defined) / OUT (unchanged)
							\%hol_inf_site_seq_of_tax_of_infile,	# key: input file name; value ; list of taxa and corresponding sequences				-> IN (defined) / OUT (unchanged)
							\$$sref_parameter_phy,					# phylip output parameter ('strict' or 'relaxed')										-> IN (defined) / OUT (unchanged)
							\%$href_hoh_info_of_infile_of_type,		# key1: infilename, key2: type of info (e.g. 'seqlength'); info type associated value	-> IN (defined) / OUT (unchanged)
							\'parsimony'							# if 'parsimony' do not set conversion counter +1										-> IN (defined) / OUT (unchanged)
				)
			}
			##############################
			
			############################################################ END PRINT OUT OF PARSIMONY INFORMATIVE SITE FOR EACH INFILE
		}
		########################################################################################## END SUBROUTINES OF &print_out
	}
	
	sub end{
		
		my $href_outfile_name	= $_[0] ;	# Outfilename of output option							-> IN (defined) / OUT (unchanged)
		my $aref_parameter_inf	= $_[1] ;	# List of parameter options of info print out			-> IN (defined) / OUT (unchanged)
		my $aref_parameter_phy	= $_[2] ;	# List of parameter options of PHYLIP print OUT			-> IN (defined) / OUT (unchanged)
		my $aref_parameter_nex	= $_[3] ;	# List of parameter options of NEXUS print OUT			-> IN (defined) / OUT (unchanged)
		my $aref_parameter_con	= $_[4] ;	# List of consensus options								-> IN (defined) / OUT (unchanged)
		my $aref_parameter_fil	= $_[5] ;	# List of of file handling								-> IN (defined) / OUT (unchanged)
		my $aref_parameter_fas	= $_[6] ;	# defined parameter option of FASTA print OUT			-> IN (defined) / OUT (unchanged)
		my $aref_parameter_3rd	= $_[7] ;	# List of third position handling						-> IN (defined) / OUT (unchanged)
		my $aref_parameter_ryc	= $_[8] ;	# List of RY coding										-> IN (defined) / OUT (unchanged)
		my $aref_parameter_tra	= $_[9] ;	# List of sequence translation options					-> IN (defined) / OUT (unchanged)
		my $aref_parameter_par	= $_[10];	# Print parsimonious sites as extra msa file			-> IN (defined) / OUT (unchanged)
		my $aref_parameter_ren	= $_[11];	# Rename taxon names of given ifiles					-> IN (defined) / OUT (unchanged)
		my $aref_parameter_prt	= $_[12];	# Print partition files for concatenated data			-> IN (defined) / OUT (unchanged)
		my $aref_parameter_mis	= $_[13];	# Replacement code of missing gene sequences 			-> IN (defined) / OUT (unchanged)
		my $sref_structure_seq	= $_[14];	# defined if structure sequence present in infiles		-> IN (defined) / OUT (unchanged)
		
		
		##############################
		## Print CLOSE SCRIPT header
		
		unless ( $aref_parameter_ren->[0] eq 'NO'	){  print "\n\n\tSEQUENCE RENAMING Process:\n\tDone"	}
		
		
		unless ( $aref_parameter_tra->[0] eq 'NO'	){
			
			 print "\n\n\tSEQUENCE Translation Process:" ;
			 
			 if ( $aref_parameter_fil->[0] =~ /Supermatrix/i	){
				
				if ( $aref_parameter_tra->[0] eq 'NUC to AA'		){ print "\n\tNUC States translated to AA States In Supermatrix Sequences"	}
				if ( $aref_parameter_tra->[0] eq 'AA to NUC'		){ print "\n\tAA States translated to NUC States In Supermatrix Sequences"	}
			}
			
			if ( $aref_parameter_fil->[0] =~ /Convert/i	){
				
				if ( $aref_parameter_tra->[0] eq 'NUC to AA'		){ print "\n\tNUC States translated to AA States In Single Infiles"	}
				if ( $aref_parameter_tra->[0] eq 'AA to NUC'		){ print "\n\tAA States translated to NUC States In Single Infiles"		}
			}
		}
		
		
		unless ( $aref_parameter_con->[0] eq 'NO'	){
			
			print "\n\n\tSEQUENCE Consensus Process:" ;
			
			if ( $aref_parameter_fil->[0] =~ /Supermatrix/i	){
				
				if ( $aref_parameter_con->[0] eq 'Freq'		){ print "\n\tDefined Consensus Blocks (Frequency) Generated Of Supermatrix Sequences"	}
				if ( $aref_parameter_con->[0] eq 'Maj'		){ print "\n\tDefined Consensus Blocks (Majority) Generated Of Supermatrix Sequences"	}
				if ( $aref_parameter_con->[0] eq 'Strict'	){ print "\n\tDefined Consensus Blocks (Strict) Generated Of Supermatrix Sequences"	}
			}
			
			if ( $aref_parameter_fil->[0] =~ /Convert/i	){
				
				if ( $aref_parameter_con->[0] eq 'Freq'		){ print "\n\tDefined Consensus Blocks (Frequency) Generated In Single Infiles"	}
				if ( $aref_parameter_con->[0] eq 'Maj'		){ print "\n\tDefined Consensus Blocks (Majority) Generated In Single Infiles"		}
				if ( $aref_parameter_con->[0] eq 'Strict'	){ print "\n\tDefined Consensus Blocks (Strict) Generated In Single Infiles"		}
			}
		}
		
		
		unless ( $aref_parameter_ryc->[0] eq 'NO'	){
			
			print "\n\n\tRY CODING Process:" ;
			
			if ( $aref_parameter_fil->[0] =~ /Supermatrix/i	){
				
				if ( $aref_parameter_ryc->[0] eq 'All'		){ print "\n\tSupermatrix Sequences Completely RY Recoded (If NUC Data)"	}
				if ( $aref_parameter_ryc->[0] eq '3rd'		){ print "\n\tSupermatrix Sequences At Each 3rd Position RY Recoded (If NUC Data)"	}
			}
			
			if ( $aref_parameter_fil->[0] =~ /Convert/i	){
				
				if ( $aref_parameter_ryc->[0] eq 'All'		){ print "\n\tSingle Infiles Completely RY Recoded (If NUC Data)"	}
				if ( $aref_parameter_ryc->[0] eq '3rd'		){ print "\n\tSingle Infiles At Each 3rd Position RY Recoded (If NUC Data)"	}
			}
		}
		
		
		if ( $aref_parameter_3rd->[0] eq 'Reject'	){
			
			print "\n\n\t3rd Position Exclusion:" ;
			
			if ( $aref_parameter_fil->[0] =~ /Supermatrix/i	){  print "\n\t3rd Positions Excluded in Supermatrix (If NUC Data)"	}
			if ( $aref_parameter_fil->[0] =~ /Convert/i		){  print "\n\t3rd Positions Excluded in Single Infiles (If NUC Data)"	}
		}
		
		
		if ( $aref_parameter_fil->[0] =~ /Supermatrix/i	){
			
			print "\n\n\tSEQUENCE Concatenation Process:" ;
			
			if ( $aref_parameter_fas->[0] eq 'YES'		){ print "\n\tConcatenated Supermatrix FILE Printed As FcC_supermatrix.fas"			}	 
			if ( $aref_parameter_phy->[0] eq 'STRICT'	){ print "\n\tConcatenated Supermatrix FILE Printed As FcC_supermatrix.phy (Strict)"	}
			if ( $aref_parameter_phy->[0] eq 'RELAXED'	){ print "\n\tConcatenated Supermatrix FILE Printed As FcC_supermatrix.phy (Relaxed)"	}
			if ( $aref_parameter_nex->[0] eq 'BLOCK'	){ print "\n\tConcatenated Supermatrix FILE Printed As FcC_supermatrix.nex (Block)"		}
			if ( $aref_parameter_nex->[0] eq 'MrBAYES'	){ print "\n\tConcatenated Supermatrix FILE Printed As FcC_supermatrix.nex (MrBayes)"	}
			if ( $aref_parameter_mis->[0] eq 'Missing'	){ print "\n\tLacking gene sequences filled with missing characters (N|X)"				}
			if ( $aref_parameter_mis->[0] eq 'Indel'	){ print "\n\tLacking gene sequences filled with indel characters (-)"					}
		}
		
		
		if ( $aref_parameter_fil->[0] =~ /Convert/i	){
			
			print "\n\n\tFILE Conversion Process:" ;
			
			if ( $aref_parameter_fas->[0] eq 'YES'		){ print "\n\tDefined INFILES Converted To *.fas FASTA" 			} 
			if ( $aref_parameter_phy->[0] eq 'STRICT'	){ print "\n\tDefined INFILES Converted To *.phy PHYLIP (Strict)"	}
			if ( $aref_parameter_phy->[0] eq 'RELAXED'	){ print "\n\tDefined INFILES Converted To *.phy PHYLIP (Relaxed)"	}
			if ( $aref_parameter_nex->[0] eq 'BLOCK'		){ print "\n\tDefined INFILES Converted To *.nex NEXUS  (Block)"	}
			if ( $aref_parameter_nex->[0] eq 'MrBAYES'	){ print "\n\tDefined INFILES Converted To *.nex NEXUS  (MrBayes)"	}
		}
		
		
		unless ( $aref_parameter_par->[0] eq 'NO'	){
			
			print "\n\n\tPrint OUT Of Parsimony Sites:" ;
			 
			if ( $aref_parameter_fil->[0] =~ /Supermatrix/i	){ print "\n\tPrint OUT Of Parsimony Informative Sites Of Supermatrix Sequences"	}
			
			if ( $aref_parameter_fil->[0] =~ /Convert/i	){ print "\n\tPrint OUT Of Parsimony Informative Sites Of Infiles Sequences"	}
		}
		
		
		unless ( $aref_parameter_prt->[0] eq 'NO'	){
			
			print "\n\n\tPartition File Print OUT:" ;
			
			if ( $aref_parameter_fil->[0] =~ /Supermatrix/i	){
				
				unless ( ( $aref_parameter_fas->[0] eq 'NO' ) && ( $aref_parameter_phy->[0] eq 'NO' ) ){ print "\n\tRAxML Partition File Printed As FcC_supermatrix_partition.txt"	}
				if ( $aref_parameter_nex->[0] eq 'BLOCK'		){ print "\n\tNEXUS (Block) Partition File Printed As FcC_supermatrix_partition.nex"	}
				if ( $aref_parameter_nex->[0] eq 'MrBAYES'	){ print "\n\tNEXUS (MrBayes) Partition File Printed As FcC_supermatrix_partition.nex"	}
			}
			
			else{
				
				print "\n\tCannot Print Partition File For Single Infile(s)"
			}
		}
		
		
		if ( $aref_parameter_inf->[0] eq 'YES'	){
			
			print "\n\n\tFILE INFO Print OUT:" ;
			
			if ( $aref_parameter_fil->[0] =~ /Supermatrix/i	){
				
				print "\n\tExtensive SEQUENCE Info Of Supermatrix Printed To ", $href_outfile_name->{info} ;
			}
			
			print "\n\tExtensive SEQUENCE Info Of INFILES Printed To ", $href_outfile_name->{info} ;
			
			
			if ( ( $aref_parameter_fil->[0] =~ /Supermatrix/i	) && ( $$sref_structure_seq ) ){ 
				
				print "\n\tExtensive Secondary STRUCTURE Info Of Supermatrix Printed To "	, $href_outfile_name->{structure} ;
				print "\n\tBasic Secondary STRUCTURE Info Of Supermatrix Printed To "		, $href_outfile_name->{info}
			}
			
			
			if ( ( $aref_parameter_fil->[0] =~ /Convert/i		) && ( $$sref_structure_seq ) ){
				
				print "\n\tExtensive Secondary STRUCTURE Info Of Infiles Printed To "	, $href_outfile_name->{structure} ;
				print "\n\tBasic Secondary STRUCTURE Info Of Infiles Printed To "		, $href_outfile_name->{info} ;
			}
		}
		
		else{ 
			
			if ( $aref_parameter_fil->[0] =~ /Supermatrix/i	){ print "\n\tBasic SEQUENCE Info Of Supermatrix Printed To "	, $href_outfile_name->{info} }
			if ( $aref_parameter_fil->[0] =~ /Conver/i		){ print "\n\tBasic SEQUENCE Info Of INFILES Printed To "		, $href_outfile_name->{info} }
		}
		
		
		print "\n";
		
		TIMER:
		# set timer
		my ( $user, $system, $cuser, $csystem ) = times ;
		
print <<TIME;
		
			***  time used: $user sec  ***
		
TIME
		
		print	"\t------------------------------------------------------------\n",
				"\t#### FASconCAT Log Off ! ####\n\n" ,
		;
		
		exit
	}
}

sub commandline{
	
	my $sref_print_text = $_[0] ;
	
	
	##############################
	## print text
	print	"\n\t", $$sref_print_text, 
			"\n\tCOMMAND:\t "  ;
	##############################
	
	
	##############################
	## READ IN and return user command
	chomp ( my $sub_answer_opening = <STDIN> );
	
	print  "\n\t------------------------------------------------------------\n" ;
	
	return $sub_answer_opening ;
	##############################
}









