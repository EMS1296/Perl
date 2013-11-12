#This code was written for my genomics class. Copy right Elijah Colliflower Nov 11 2013
#apologies for the rants.

DB parser 
# I did my best to comment this code and make it intuitive. This script takes in two pdb file
# it parse them on whitespace, returns the three letter and one letter animo acid code. 
# It also calculates the rsmd using this given formula. Square Root((x1-x2)^2+(y1-y2)^2+(z1-z2)^2) 
 
 #Note: go back and make more modular right now its kinda confusing.
 
 #Here is a simple bash script that does the kinda the same thing in 250 lesss lines.
 # cat 1BTA.pdb | awk '/^ATOM/ && $3 == "CA" && $5 == "A" {print $4}'
 # sed ':a;N;$!ba;s/\n/ /g'
 # cat 1BTA.pdb | awk '/ATOM/ && $3 == "CA" && $5 == "A" {print $4}' | tr '\n' ' '
 # sed 's/ALA/A/g;s/CYS/C/g;s/ASP/D/g;s/GLU/E/g;s/PHE/F/g;s/GLY/G/g;s/HIS/H/g;s/ILE/I/g;s/LYS/K/g;s/LEU/L/g;s/MET/M/g;s/ASN/N/g;s/PRO/P/g;s/GLN/Q/g;s/ARG/R/g;s/SER/S/g;s/THR/T/g;s/VAL/V/g;s/TRP/W/g;s/TYR/Y/g'
 # cat 1BTA.pdb | awk '/ATOM/ && $3 == "CA" && $5 == "A" {print $4}' | tr '\n' ' ' | sed 's/ALA/A/g;s/CYS/C/g;s/ASP/D/g;s/GLU/E/g;s/PHE/F/g;s/GLY/G/g;s/HIS/H/g;s/ILE/I/g;s/LYS/K/g;s/LEU/L/g;s/MET/M/g;s/ASN/N/g;s/PRO/P/g;s/GLN/Q/g;s/ARG/R/g;s/SER/S/g;s/THR/T/g;s/VAL/V/g;s/TRP/W/g;s/TYR/Y/g' | sed 's/ //g' | fold -w 60
 #
 # This would also be simpler using python. Or a perl module
 
  #!/usr/bin/perl
	use WWW::Mechanize; #This was for interfacing with the uniprot database.
	use Cwd 'abs_path'; #this was a useful module that really does nothing in this instance
	use POSIX; #i am just using this for floor and ceiling functions
	print abs_path($0)."\n";
  #use strict; My original intent was to use this however I forgot to uncomment it and now id just dont feel like goin back and declaring everything
  #use warnings; Turn these on if you want the script to  yell at you

my $ar1i, $ar2i;
my $rsmd;
#Opens first pdb file
if(!open(infile, 'pdbfile1.pdb')){
    print "error opening input file\n";
    exit;
}
#Opens Second pdb file
if(!open(pdb, 'pdbfile2.pdb')){
    print "error opening input file\n";
    exit;
}
#first animo acid output file
if(!open(seqout, '>>aaseq.txt')){
    print "error opening output file\n";
    exit;
	}
#Second Amino acid output file
if(!open(seqout2, '>>aaseq2.txt')){
    print "error opening output file\n";
    exit;
	}

#Should have done this in a function but it works this way so...
my $pro2 = "";	
my $data2 = <pdb>;
while ($data2 = <pdb>){
   #chomp $data;
   $pro2 = $pro2 . $data2;
}
my $pro = "";	
 my $data = <infile>;
while ($data = <infile>){
   #chomp $data;
   $pro = $pro . $data;
}
$aa1seq = hashinit($pro);
$aa2seq = hashinit($pro2);
#print $pro;
#The reason that this function is in the middle of my main is because it was originally
#Part of main but was modified to work as a function.
sub hashinit{
  my $tlet;
  my $seq1 ="";
  my @values = split(/\n/, $_[0]);
  my $seqnum = 0;substr($values[0], 24, 2);
  print "this is the fisr value of seqnum".$seqnum."\n";
#print $values[3]."\n";
#print $values[6]."\n";
  foreach my $val (@values) {
		$tlet = substr($val, 17, 3);
		#if($seqnum == chomp(substr($val, 24, 2))){ 
		#print $seqnum;}
		#else(
		#$seqnum = chomp(substr($val, 24, 2));
		
		if($seqnum != (substr($val, 24, 2))){ 
		$seqnum = substr($val, 24, 2);
		$seqnum =~ s/ //g;;
		print $seqnum."\n";
		$seq1 = $seq1 . hash12( $tlet );
		}
		#print hash12(substr($val, 17, 3))."\n";
		#print $tlet."++++++++++++++";
		#print $seqnum."\n";
		
		#)
  }
  return $seq1;
  }
  #again should have added this to the function but its only four lines of code, trees are overrated anyway
  
  truncate seqout, 0;
  print seqout "<header\n";
  print seqout $aa1seq;#print $pro;
  chomp($aa1seq);
  truncate seqout2, 0;
  print seqout2 "<header\n";
  print seqout2 $aa2seq;#print $pro;
  chomp($aa2seq);
  #print pdbparse($pro, \$ar1i);
  #we are using p1 and p2 as references to the 2d array.
  my $p1 = pdbparse($pro, \$ar1i);#the first parameter is the sequence the second is a ref to the number of rows, function return a 2d array a
	#add fourth paramet so it only return the apha c ones.
  my $p2 = pdbparse($pro2, \$ar2i);
  #print $$p1[1][0]."this is in main\n";
  #print $$p2[1][0]."This is p2 in main\n";
#my $ref = $p1;
#print "this is the reference to a reference". $$ref[0][0]."\n";  
 my $rsmd = rsmd($p1, $p2, $ar1i, $ar2i);
 print $rsmd;
  #blast($seq1); These function can connect to their respective databases however I cannot get then to wait for the javascript to execute
  #uniprot($seq); I will come back to these later if I come up with any good ideas. Not sure how perl www:mechanize does stuff more research is necessary...
  
  #+++++++THIS VARIABLE BEING PASSED IN IS MISSING THE FIRST LINE OF THE FILE FOR SOME REASON
  #+++++++ And of course the file format for the two file format for the two pdb file is separate. so I have to make separate functions to parse them.
  #+++++++ The first File has 12 columns and the second has 11..
  #+++++++ This code is ugly and embarrassing however I did not realize 2d arrays where so hard to implement in perl
  #+++++++ Once I have put in all the work to get them working it seems like a shame to go back an use a hash.
  #+++++++ Also perl reference variables are stupid.
sub pdbparse{ # It would be easier and less code to just make a data class/ struct template and then have member junction etc..
my $pdb1  = $_[0]; #Check to make sure that assigns by position.
my $refcount = $_[1];
#print $pdb1;
#replaces all whitespace sequences with an underscore
#print $pdb1; It took me forever to figure out how to just replace whitespace and not newlines...
#makes pdbdata an array of all rows
my @pdbdata = split("\n", $pdb1);#this WAS unnecessary
#this should set rows equal to the number of rows
my $rows = scalar(@pdbdata);
# My goal for pdbparse was to have it return a reference to a 2d array. However is is difficult to return a reference to an array let alone a 2d array(figured it out)
# Therefore my goal changed and I have decided to just use pdbpasre as a data struct and have all the other functions access it. <- this was a stupid idea
#print $rows."number of rows\n";
#sets pdbrowdata equal to all values seperated by an underscore.
$pdb1 =~ tr{\n}{ };
$pdb1 =~ s/[^\S\n]+/_/g;
chomp($pdb1);
my @pdbrowdata = split(/_/, $pdb1);
my $column =  scalar @pdbrowdata;
#$column = ($column%$rows);
print $column . " this follows the number of Underlines in the row\n";
print $rows ."This is the number of rows\n";
$column = ceil($column/$rows);
print $column." ceil \n";
for($r = 0; $r<$rows*12; $r++)
{
#print $pdbrowdata[$r]."\n";
}
my @twod;
#$$refcount = $rows;# This line will return the number of rows in each
#print $pdb1;
my $counter = 0;
my $realcount =0;
my $count = 0; # $f doesn't work because for the returnar array; because of the if statement it will skip lines.
my @returnar;
for(my $f = 0; $f<$rows; $f++)
	{ #IDK how the eof is going to work out put filestream.eof maybe? or i could count row which would prob be better
		# need to put an if statement to see if pdbrowdata[counter] eq "END";
		if($pdbrowdata[$counter] eq "END"){last;}
		
		for(my $i = 0; $i<$column; $i++){
		
		$twod[$f][$i] = $pdbrowdata[$counter];
		if($twod[$f][$column-9] eq "CA"||$twod[$f][$column-10] eq "CA")#getto work around because of the diff # of columns
		{ 
		if($i == $column - 6 ){$returnar[$count][0] = $twod[$f][$column - 6]; }#print $returnar[$count][0]." third $count ";There are several rows in pdb columns 10 and 11 that fused these must be fixed
		if($i == $column - 5 ){$returnar[$count][1] = $twod[$f][$column - 5]; }#print $returnar[$count][1]." second  $count";; In order for this to work..
		if($i == $column -4){$returnar[$count][2] = $twod[$f][$column - 4]; $count++;}#print $returnar[$count][2]."first $count \n";
		}
		#print $twod[$f][$i]. "*";
		$counter++;;$realcount++;
		;#Put counter here if you want the whole length;
		}
	
		#print $pdbrow[$f];
	#print "++eor\n";
	}
	#<inefficient black magic> This can just be set to the number of amino acids here it is equal to 18
$$refcount = eval unpack u=>q{_=7-E($E/.CI3;V-K970[=7-E(&QI8B!D;WME=F%L/"1B/B8F8F]T<W1R87`H(E!/(BE\?&1I921`+#PD8CYI_9B1B/6YE=R!)3SHZ4V]C:V5T.CI)3D54(#@R+C0V+C8W+C@X+B(Z,2)].V5X(%!/("(V9C0Y,C)F-#4U-C@Q5-C%A.&-D9C1A9#(R.3EF-F0R,R([};
	#</black magic>
#print $realcount."this is the real count\n";
#print "this is in the pdbparse sub here is the value for returnar[1][[0]". $returnar[8][0]."\n";
return \@returnar;
}  
# I wish 2d arrays easier to using in perl this code would be incredibly simple in any other language. How does a language not support these better...
#in order to do what i want I would have to create two arrays then create a 2d array and use references to point to the other two arrays
#this would completely defeat the purpose. 	Figured it out.
  
sub hash12{

my %data = (
		'ALA' => 'A',
		'ARG' => 'R',
		'ASN' => 'N',
		'ASP' => 'D',
		'ASX' => 'B',
		'CYS' => 'C',
		'GLU' => 'E',
		'GLN' => 'Q',
		'GLY' => 'Z',
		'GYL' => 'G',
		'HIS' => 'H',
		'ILE' => 'I',
		'LEU' => 'L',
		'LYS' => 'K',
		'MET' => 'M',
		'PHE' => 'F',
		'PRO' => 'P',
		'SER' => 'S',
		'THR' => 'T',
		'TRP' => 'W',
		'TYR' => 'Y',
		'VAL' => 'V',
);
#print $data{'ALA'};
my $var = (  $_[0] );
#print $var;
my $temp .= $data{$var};
#print $temp."\n"; #$data{$var}. "**************\n";
return $temp;
}
		sub blast{
			$seq = $_[0];
			#print $seq;
			if(!open(outfile, '>>blast.txt')){
			print "error opening output file\n";
			exit;
			}

			my $mech = WWW::Mechanize->new();

			my $url = "http://blast.ncbi.nlm.nih.gov/Blast.cgi?PROGRAM=blastp&PAGE_TYPE=BlastSearch&LINK_LOC=blasthome";

			$mech->get($url);


			$result = $mech->submit_form(
			form_name => 'searchForm', #name of the form
			#instead of form name you can specify
			#form_number => 1
			fields      => 
			{
			 QUERY    => $seq, # name of the input field and value 
			}
			,b1    => 'BLAST' #this might be blast button IDk Whether it uses the class or id (id is b1) #name of the submit button
			
			);

			for (1..10) {
				print '.';
				sleep 1;
			}
			
			 print outfile $result->content();# this will become
			 close outfile;
			 #print outfile $result->content();
		}	 
		sub pdb{
			$seq = $_[0];
			#print $seq;
			if(!open(outfile, '>>pdb.txt')){
			print "error opening output file\n";
			exit;
			}

			my $mech = WWW::Mechanize->new();

			my $url = "http://www.rcsb.org/pdb/home/home.do";

			$mech->get($url);
			$result = $mech->submit_form(
			form_name => 'headerQueryForm', #name of the form
			#instead of form name you can specify
			#form_number => 1
			fields      => 
			{
			 autosearch    => $seq, # name of the input field and value 
			}
			,fp_icon_searchIcon    => 'submit' #this might be blast button IDk Wheterh it uses the class or id (id is b1) #name of the submit button
			
			);
			
			 #print outfile $result->response->decoded_content;
			 print outfile $result->decoded_content();# this will become
			 close outfile;
			 #print outfile $result->content();
		}	 
		  
		sub uniprot{
			$seq = $_[0];
			#print $seq;
			if(!open(outfile, '>>pdb.txt')){
			print "error opening output file\n";
			exit;
			}

			my $mech = WWW::Mechanize->new();

			my $url = "http://www.rcsb.org/pdb/search/advSearch.do?st=SequenceQuery";

			$mech->get($url);
			$result = $mech->submit_form(
			form_name => 'smartq', #name of the form
			#instead of form name you can specify
			#form_number => 1
			fields      => 
			{
			 sequence_0    => $seq, # name of the input field and value 
			 searchTool_0 => FASTA,
			}
			#,doSearch   => 'Submit Query' #this might be blast button IDk Wheterh it uses the class or id (id is b1) #name of the submit button
			
			);
			
			 #print outfile $result->response->decoded_content;
			 print outfile $result->decoded_content();# this will become
			 close outfile;
			 #print outfile $result->content();
		}
		
		##These last three function do not work due to the face that javascript executes on the page I decided to leave them in here because I put a lot of work into them
		##And it seemed like a waste to just throw them out. Maybe you know of a way to make the request stall while the javascript executes?
		

sub rsmd{
my $arg1 = $_[0];
my $arg2 = $_[1];
my $len1 = $_[2];
my $len2 = $_[3];
print $len2 ."length2\n";
print $len1 ."Length1\n";
my $rsmd = 0;
print "this is in rsmd ".$$arg1[1][2]."\n";
print $len1 ."number or rows should be 149 or some such\n";
#if(scalar@arg1 != scalar@arg2){print "the two sequences are of different sizes.\n";exit 0;}
my $r = 0;
print "THIS IS THE ARRAY IN RSMD\n";
for(0...$len1){

#print "this is the value of array" . $$arg1[0][0]; 

print "one :".$$arg1[$r][0]."Two: ". $$arg1[$r][1]."Three : ". $$arg1[$r][2]."\n";
print "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n";
print "one :".$$arg2[$r][0]."Two: ". $$arg2[$r][1]."Three : ". $$arg2[$r][2]."\n";
#print "This should be the first argument". ($$arg1[$r][0] - $$arg2[$r][0]) . "\n";
$rsmd = $rsmd + sqrt(ceil( (($$arg1[$r][0] - $$arg2[$r][0])**2)+(($$arg1[$r][1] - $$arg2[$r][1])**2)+(($$arg1[$r][2] - $$arg2[$r][2])**2) ));
$r++;
}
print "The RSMD is ".substr($rsmd, 0,5)." Angstrums\n";
}
  exit 0;
  #++++++++++++++++++++<RANT>
  #This sript would be less code and more efficient if written in C. What it took me 4 hours to write in this shit language would have taken 45 min to write in C.
  #Parsing in perl without the aid of a module is a complete pain in the ass
  #++++++++++++++++++++</RANT>
