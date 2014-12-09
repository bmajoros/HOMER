#!/usr/bin/perl
#**************************************************************
# auto-trainer.pl
#
# Distributed as part of the "Gene Finding for Dummies" book.
#
#**************************************************************
use strict;
$|=1;

my $maxTrain=1000; ### not sure what to do about this...

#**************************************************************
# Make sure this is the correct directory for training
#**************************************************************
my $cwd=`pwd`;
print "
Hello!  I am the auto-trainer for the Unveil gene finder.
I will be assisting you to train Unveil for a specific organism.

You should be running this script in the same directory where
the training files are located (which should NOT be the same
directory where Unveil is installed).

The current directory is:

$cwd

Is this the directory containing the training files for the
desired organism?

";

my $input=YorN();
if($input eq "N")
  {
    print "
You have indicated that this is not the correct directory.
Please change into the correct directory and re-run this
script from there.

";
    exit;
  }


#**************************************************************
# Set environment variables
#**************************************************************
print "
Now we need to set some environment variables.  Please type the
full path of the directory where Unveil is installed:

";
my $dir=<STDIN>;chop $dir;
if($dir eq $cwd) 
  {
    die "
Uh-oh!  It appears that the training directory is the same as
the unveil installation directory.  This is not allowed.  Please
move the training files into a different directory, cd into that
directory, and re-run this script from there.

";
  }
print "
OK, I'm going to set the PERLLIB, PERL5LIB, and PATH variables now...
";
$ENV{"PERLLIB"}=$ENV{"PERL5LIB"}="$dir/perl:$dir/perlib";
$ENV{"PATH"}.=":$dir:$dir/perl";
print "
Set.  I'm also going to write these into a file so you know how to
set these variables yourself later.  I'll put the definitions in
a file called \"env.txt\"...
";
while(!open(OUT,">env.txt"))
  {
    print "
Uh-oh.  I wasn't able to open that file in write-mode.  You may need
to change the permissions of the current directory.  Try something
like: \"chmod a+wx .\".  When you think you've fixed the problem I'll
try again.

Press ENTER when you're ready for me to try again...";
<STDIN>;
  }
print OUT "
csh
setenv PERLLIB $dir/perl:$dir/perlib
setenv PERL5LIB $dir/perl:$dir/perlib
setenv PATH \$\{PATH\}:$dir:$dir/perl
";
close(OUT);
print "
Done.  Now let's create a symbolic link to the unveil directory...
";

#**************************************************************
# Create symbolic link to unveil directory
#**************************************************************
while(-e "unveil")
  {
    my $type=`file unveil`;
    if($type!~/symbolic\s*link/ || !-e "unveil/model-combiner")
      {
	die "
Uh-oh.  I need to create a symbolic link from the current directory
to the Unveil installation directory.  That symbolic link needs to
be named \"unveil\".  It seems there is already a file called
\"unveil\" in the current directory.  Please remove or rename that
file and let me try again.

Press ENTER when you are ready for me try again...";
	<STDIN>;
      }
    else {last}
  }
system("ln -s $dir unveil") unless(-e "unveil");
while(!-e "unveil/model-combiner")
  {
    die "
Uh-oh.  I tried to create a symbolic link called \"unveil\" to the
Unveil installation directory.  For some reason, this operation 
failed.  Please check that you have write permission in the
current directory and allow me to try again.

Press ENTER when you are ready for me to try again...";
    <STDIN>;
    system("ln -s $dir unveil");
  }

#**************************************************************
# Extract the training features
#**************************************************************
print "
OK, now I'm ready to process the training files.  You need to
provide two files: (1) a GFF file containing exon coordinates,
and (2) a FASTA file to which those coordinates refer.  Both of
those files must be in the current directory.

What is the name of the GFF file containing the training data?
You don't need to specify the full path -- just the filename.

==> ";
my $trainGff=<STDIN>;chop $trainGff;
while(! -e $trainGff)
  {
    print "
Uh-oh.  I can't find that file in the current directory.  Are
you sure it's there?  Have you typed the name correctly?  Please
verify that the file is present and I will try again.

Name of GFF file ==> ";
    $trainGff=<STDIN>;chop $trainGff;
  }

print "
OK, now give me the name of the associated FASTA file:

==> ";
my $trainFasta=<STDIN>;chop $trainFasta;
while(! -e $trainFasta)
  {
    print "
Uh-oh.  I can't find that file in the current directory.  Are
you sure it's there?  Have you typed the name correctly?  Please
verify that the file is present and I will try again.

Name of FASTA file ==> ";
    $trainFasta=<STDIN>;chop $trainFasta;
  }
system("unveil/get-training-files.pl $trainGff $trainFasta $maxTrain");

#**************************************************************
# Build HMM templates
#**************************************************************
print "
Splendid.  Now I will attempt to build the default HMM template
files...
";
system("unveil/build-hmm-templates.pl");
while(!-e "exon.hmms" || !-e "submodels.txt")
{
  print "
Uh-oh.  I ran the unveil/build-hmm-templates.pl script, but it
did not create all of the files I was expecting.  Did you see
any error messages when I ran that script?  If so, please try
to address the problem and I will try to run the script again.

Press ENTER when you are ready for me to try again...";
  <STDIN>;
}

#**************************************************************
# Train submodels
#**************************************************************
print "
OK, now I am going to try to train the submodels of the HMM.
This may take a long time, depending on the speed of your
computer.  You may want to observe the output of this process
to see if any errors occur.  You may ignore warnings about
decreases in log-likelihood (they are due to small rounding
errors in the computer).

Press ENTER when you are ready to begin training...";
<STDIN>;
system("unveil/train-submodels.pl");
unless(-e "acceptor-frame-shift.hmm" && -e "downstream-intergenic.hmm"
       && -e "stop-codon.hmm" && -e "acceptor.hmm" && -e "exon.hmm"
       && -e "donor-frame-shift.hmm" && -e "intron.hmm" 
       && -e "upstream-intergenic.hmm" && -e "donor.hmm" 
       && -e "start-codon.hmm")
{
  die "
Uh-oh.  One or more of the *.hmm files that the training process
should have created is missing.  You may want to scroll up this
window to see if any errors occurred.  If you are able to mend
the problem, you can resume the training process by running this
script again from scratch.  Sorry.
";
}

#**************************************************************
# Combine models
#**************************************************************
print "
Great, it looks like the individual submodels have been
successfully trained.  Now we just need to combine the submodels
into a single HMM.  This should only take a second.

Press ENTER when you are ready...";
<STDIN>;
system("unveil/model-combiner metamodel.hmms submodels.txt unveil.hmm");
print "
Done.  Did you see any error messages?

";
my $yn=YorN();
while($yn eq "Y")
{
  print "
Oops.  Sorry about that.  If you think you can fix the problem,
go ahead and I will re-run the model-combiner when you are ready.

Press ENTER when you are ready...";
<STDIN>;
print "
Did you see any error messages that time?

";
my $yn=YorN();
}

#**************************************************************
# Make test set
#**************************************************************
print "
Congratulations!  The gene finder has been completely trained.
Now let's run it on the training set just as a sanity check, to
make sure nothing untoward happened during training.  Because
we're testing on the training set, the accuracy measurements
will be artificially inflated, but they will serve our purpose
for the time being (just don't cite the numbers in a publication).

I am going to begin by creating a \"chunks\" directory and
populating it with 100 test genes.

Press ENTER when you are ready...";
<STDIN>;
system("unveil/make-test-set.pl $trainGff $trainFasta 100 500");
while(!-e "chunks/1.gff" || !-e "chunks/1.fasta")
{
  if(!-e "chunks")
  {
    print "
Uh-oh.  I wasn't able to create a \"chunks\" subdirectory.  Perhaps
you need to set the permissions using \"chmod a+wx .\"?  I'll
try again when you have fixed the problem.

Press ENTER when you are ready...";
  }
  else
  {
    print "
Uh-oh.  That didn't seem to work.  The \"chunks\" directory appears
to have been created successfully, but there don't appear to be
any chunks in there.  Perhaps you need to change the permissions of
that directory using \"chmod\"?  I'll try again when you think you 
have fixed the problem.

Press ENTER when you are ready...";
  }
  <STDIN>;
  system("unveil/make-test-set.pl $trainGff $trainFasta 100 500");
}

#**************************************************************
# Run unveil on the test set
#**************************************************************
print "
Divine.  Now we just have to run the gene finder on the chunks
and evaluate the results.  First, to run the gene finder.

Press ENTER when you are ready...";
<STDIN>;
system("unveil/run.pl unveil.hmm submodels.txt");
while(!-e "out/1.gff")
{
  if(!-e "out")
  {
    print "
Uh-oh.  I wasn't able to create an \"out\" subdirectory for the
output files.  Perhaps you need to set the permissions using 
\"chmod a+wx .\"?  I'll try again when you have fixed the problem.

Press ENTER when you are ready...";
  }
  else
  {
    print "
Uh-oh.  That didn't seem to work.  The \"out\" directory appears
to have been created successfully, but there don't appear to be
any output files in there.  Perhaps you need to change the permissions
of that directory using \"chmod\"?  I'll try again when you think 
you have fixed the problem.

Press ENTER when you are ready...";
  }
  <STDIN>;
  system("unveil/run.pl unveil.hmm submodels.txt");
}

#**************************************************************
# Evaluate the test predictions
#**************************************************************
print "
Excellent.  Now we can evaluate the results:

";
system("unveil/evaluate.pl");
print "

You should see some accuracy measurements above.  If you do not,
then something has gone amiss.  If you are able to correct the
problem, you can re-run the evaluation script yourself using the
command: unveil/evaluate.pl

Finally, remember that you must execute the commands in \"env.txt\"
before running the gene finder, to set the environment variables.

Have a nice day!
";

#=============================================================
sub YorN
  {
    print "(Y/N) ==> ";
    while(1)
      {
	my $reply=<STDIN>;
	$reply=~s/\s+//g;
	$reply="\U$reply";
	if($reply eq "Y" || $reply eq "N") {return $reply}
	print "I didn't understand that response.  Please type Y or N: ";
      }
  }
#=============================================================
#=============================================================
#=============================================================
#=============================================================
#=============================================================
#=============================================================


