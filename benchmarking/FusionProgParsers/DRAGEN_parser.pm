package DRAGEN_parser;

# see Illumina support Dragen doc

use strict;
use warnings;
use Carp;


=DRAGEN_format

0       FusionGene
1       Score
2       LeftBreakpoint
3       RightBreakpoint
4       Filter
5       SplitScore
6       NumSplitReads
7       NumSoftClippedReads
8       NumSoftClipReadsGene1
9       NumSoftClippedReadsGene2,
10      NumPairedReads,
11      NumRefSplitReadsGene1,
12      NumRefPairedReadsGene1,
13      NumRefSplitReadsGene2,
14      NumRefPairedReadsGene2,
15      AltToRef,
16      UniqueAlignmentsGene1,
17      UniqueAlignmentsGene2,
18      MaxMapqGene1,
19      MaxMapqGene2,
20      CoverageBasesGene1,
21      CoverageBasesGene2,
22      DeltaExonBoundaryGene1,
23      DeltaExonBoundaryGene2,
24      IsRestrictedGene1,
25      IsRestrictedGene2,
26      IsEnrichedGene1,
27      IsEnrichedGene2,
28      CisDistance,
29      BreakpointDistance,
30      GenePairHomologyEval,
31      AnchorLength1,
32      AnchorLength2,
33      FusionLengthGene1,
34      FusionLengthGene2,
35      NonFusionLengthGene1,
36      NonFusionLengthGene2,
37      AdditionalGenes1,
38      AdditionalGenes2,
39      Gene1Id,
40      Gene2Id,
41      Gene1Location,
42      Gene2Location,
43      Gene1Sense,
44      Gene2Sense



0          ACAP1--CRYBB3
1          1
2          chr17:7342067:+
3          chr22:25202674:+
4          PASS
5          1761
6          1304
7          457
8          241
9          216
10         699
11         103
12         36
13         1022
14         541
15         1.72309
16         417
17         439
18         133
19         133
20         884
21         1000
22         0
23         0
24         1
25         1
26         0
27         0
28         100000000
29         100000000
30         1
31         90
32         101
33         5539
34         4686
35         9411
36         2817
37
38
39         ENSG00000072818.12
40         ENSG00000100053.10
41         IntactExon
42         IntactExon
43         True
44         True


=cut


sub parse_fusion_result_file {
    my ($file) = @_;

    my @fusions;
    
    open (my $fh, $file) or die "Error, cannot open file $file";
    my $header = <$fh>;
    while (<$fh>) {
        chomp;

        my @x = split(",");
        my ($fusion_gene_A, $fusion_gene_B) = split(/--/, $x[0]);
        if ($fusion_gene_A eq $fusion_gene_B) { next; } # no self-fusions

        my $chr_coords_A = $x[2];
        my $chr_coords_B = $x[3];
        
        my ($chrA, $coordA, $orientA) = split(/:/, $chr_coords_A);
        my ($chrB, $coordB, $orientB) = split(/:/, $chr_coords_B);

        my $junction_reads = $x[6] + $x[7];
        my $spanning_reads = $x[10];
        
        my $struct = {
            geneA => $fusion_gene_A,
            chrA => $chrA || ".",
            coordA => $coordA || ".",

            geneB => $fusion_gene_B,
            chrB => $chrB || ".",
            coordB => $coordB || ".",

            span_reads => $spanning_reads,
            junc_reads => $junction_reads,
        };
        
        push (@fusions, $struct);

    }

    close $fh;

    return(@fusions);
}


1; #EOM

