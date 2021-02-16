#!/usr/bin/env perl

use strict;
use warnings;

my $usage = "usage: $0 consolidated_edgren_predictions.dat\n\n";

my $input_file = $ARGV[0] or die $usage;


main: {

    my @progs_to_count = `cat ../progs_select.txt`;
    chomp @progs_to_count;
    
    my %progs_to_count = map { + $_ => 1 } @progs_to_count;
    
    
    open(my $fh, $input_file) or die $!;

    my %fusion_to_prog;
    my %orig_fusion_call;
    my %prognames;

    my $header = <$fh>;
    while(<$fh>) {
        chomp;
        my @x = split(/\t/);

        unless ($x[1] eq "TP") { next; }

        my $sample_name = $x[2];
        my $progname = $x[3];
        my $fusion_name = $x[6];

        my @y = split(/\|/, $fusion_name);
        
        my ($left_entry, $right_entry) = split(/--/, $y[1]);
                
        $prognames{$progname}++;
        
        my $alt_fusion_name = "$sample_name|$left_entry--$right_entry";
        $orig_fusion_call{$alt_fusion_name}->{$fusion_name} += 2; # give additional weight to original orientation.
        $fusion_to_prog{$alt_fusion_name}->{$progname}++;
        
        $alt_fusion_name = "$sample_name|$right_entry--$left_entry";
        $orig_fusion_call{$alt_fusion_name}->{$fusion_name}++;
        $fusion_to_prog{$alt_fusion_name}->{$progname}++;
                
        
    }
        


    # only use one orientation.
    my %seen;
    my %fusions_meet_min_prog_count;
    foreach my $fusion_name (keys %fusion_to_prog) {
        
        if ($seen{$fusion_name}) { next; }

        my $orig_fusion_names_href = $orig_fusion_call{$fusion_name};
        my @orig_fusion_cand_names = sort {$orig_fusion_names_href->{$b}<=>$orig_fusion_names_href->{$a}} keys %$orig_fusion_names_href;

        foreach my $orig_fusion_cand_name (@orig_fusion_cand_names) {
            $seen{$orig_fusion_cand_name} = 1;
        }
        
        my $orig_fusion_name = $orig_fusion_cand_names[0];
                
        #print "$fusion_name\t$orig_fusion_name\t$prog_count\n";

        $fusions_meet_min_prog_count{$orig_fusion_name} = 1;
        
    }
        

    ## generate report
    my @prognames = sort keys %prognames;

    print "\t" . join("\t", @prognames) . "\n";

    my @final_fusions = sort keys %fusions_meet_min_prog_count; 
    
    foreach my $fusion (@final_fusions) {
        
        my @vals = ($fusion);
        foreach my $progname (@prognames) {
            my $found = (exists $fusion_to_prog{$fusion}->{$progname}) ? 1 : 0;
            push (@vals, $found);
        }

        print join("\t", @vals) . "\n";
    }

    exit(0);
    
}
