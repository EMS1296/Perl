#!/usr/local/bin/perl
use Tk;

my $mw = new MainWindow; # Main Window

my $lab = $mw -> Label(-text=>"Enter DNA Sequence") -> pack();
my $ent = $mw -> Entry() -> pack();

my $lab = $mw -> Label(-text=>"Enter Secont DNA Sequence") -> pack();
my $ent = $mw -> Entry() -> pack();

my $but = $mw -> Button(-text => "Global Alignment", 
		-command =>\&push_button);
$but -> pack();

my $cb = $mv->Checkbutton( -text => "Show Matrix", option => value[ , . . . ] )->pack(-side => 'top');

my $textWidget=$mv->Text(
          -height=>'20',
          -width=>'10'
          )->pack();

my $but = $mw -> Button(-text => "Semi-Global Alignment", 
		-command =>\&push_button);
$but -> pack();


my $but = $mw -> Button(-text => "Local Alignment", 
		-command =>\&push_button);
$but -> pack();

MainLoop;

#This is executed when the button is pressed
sub push_button {
	$ent -> insert(0,"Hello, ");
}
