# De Bruijn cycle generator

Stimulus counter-balance is important for many experimental designs. This command-line software creates pseudo-random sequences with arbitrary levels of counter-balance. "Path-guided" de Bruijn cycles may also be created. These sequences encode a hypothesized neural modulation at specified temporal frequencies, and have enhanced detection power for BOLD fMRI experiments. These sequences are particularly well-suited for carry-over fMRI experiments.

This paper describes the use of de Bruijn cycles in neuroscience experiments in general, and fMRI adaptation experiments in particular:

  * GK Aguirre, MG Mattar, L Magis-Weinberg. (2011) [de Bruijn cycles for neural decoding](http://www.ncbi.nlm.nih.gov/pubmed/21315160). _NeuroImage_ 56: 1293-1300

The sections below describe the use of our open source software for the creation of these sequences.

## Included code

The software contains Hamiltonian Cycle code that was developed by Basil Vandegriend and Joseph Culberson at the University of Alberta (Copyright 1998). Please see their [Conditions of Use](http://webdocs.cs.ualberta.ca/~joe/Theses/HCarchive/conditions.html) page regarding the further use and dissemination of this code. In a personal communication (July 4, 2013), Dr. Vandegriend has provided authorization to release the code under the [BSD-3](http://en.wikipedia.org/wiki/BSD_licenses#3-clause_license_.28.22Revised_BSD_License.22.2C_.22New_BSD_License.22.2C_or_.22Modified_BSD_License.22.29) license.

## Help text

Calling with either ''-h'' or with no parameters produces this help text:

```
<code - debruijn -h>
Usage: debruijn [-t | -v] <k> <n>
	[<B> <guide function> <neural model>]
	[<neural model #2>] [<neural model #3>] [-eval <SOA>]
Generate random and "path guided" de Bruijn cycles
  
Default output has one line of labels, separated by commas
	-t: terse output (no delimiters)
	-v: verbose output in "necklace" format
  
  	k: number of letters (maximum is 36)
	n: length of the word (i.e., level of counterbalance)

Specify parameters for path-guided cycle
	B: number of bins for available transitions (integer between 1 and k^2)
	guide function: can be specified in three different ways -
		HRF: guide with power spectrum of the BOLD HRF  -OR-
		path to .txt file containing guide function -OR-
		[Tmin,Tmax]: range of the periods of sinusoids
		(in units of labels) to use as a guide function
	neural model: path to .txt file containing neural model for optimizing
	neural model(#2,#3): optional; detection power only reported
  	-eval: Evaluate detection power
	SOA: stimulus-onset asynchrony (in ms)

Code by: Marcelo Mattar (mattar@sas.upenn.edu)
	Dongbo Hu (dongbo@mail.med.upenn.edu)
	Hamiltonian circuit code from Basil Vandergriend
	(http://webdocs.cs.ualberta.ca/~joe/Theses/vandegriend.html)
```

## Sequences and Cycles
Full counterbalance is present with the de Bruijn _cycle_; for a linear sequence to exhibit complete counterbalance, the last n elements of the sequence must be presented prior to the start of the sequence to establish the context for the first stimulus. Moreover, it is optimal in BOLD fMRI studies to present 8-15 seconds of stimuli from the end of the sequence at the start, to allow the hemodynamic response to reach steady-state. The neural data corresponding to these initial stimuli may be dropped, resulting in a final fMRI signal that represent the entire sequence of stimuli, with the delayed and dispersed BOLD fMRI response to the neural sequence “wrapping around” from end to start. A similar "repeat and clip" technique may be used at the boundaries of scans if a lengthy de Bruijn sequence must be broken into portions for subject comfort.

## Expanded usage notes

### Output format
The de Bruijn cycle will be given as a sequence of labels. The characters ''0-9'' are used for the first 10 labels; the characters ''A-Z'' are used thereafter. The maximum *k* permitted is 36.

### Number of bins
Path-guiding selects amongst available nodes whose transition match that called for by the guide function. To provide stochasticity, the available nodes are cast into bins, and the next node in the sequence is selected from within the best matching bin. A large number of bins produces a path-guided sequence of greater initial fidelity, but with decreasing precision as the sequence progresses. Conversely, a small number of bins produces a consistently less precise approximation of the guide function.

Better performance is provided when the number of nodes in each bin is roughly equivalent. For this reason, we recommend selecting B such that [k<sup>2</sup> modulo B] = 0.

### Neural Model
The *neural matrix* represents the expected neural response to stimuli (or transitions between stimuli). Using this information, "path guided" sequences may be created which are optimized for the detection of the hypothesized neural effect.

The parameter provides the path to a *k* x *k* matrix of positive floating point values in a plain text file. Columns of the matrix are delimited by spaces, and rows by carriage returns. Example neural matrices are included in the repository in .txt files.

**NOTE** The code currently requires that the neural matrix file end with the last row of the matrix. Ending the file with an additional newline or carriage return character causes the program to fail with the message: ''The neural model file must have a nxn matrix of floating points''.

For some stimuli (or transitions), the entry in the neural matrix may be left undefined. For example, from a target stimulus or for the stimulus following a "null" trial. For these transitions, a value of ''-1'' may be entered in the neural model matrix. These undefined transitions are randomly distributed throughout the sequence.

### Direct Effects
The *direct effect* is the neural response to a stimulus itself (as opposed to stimulus context). Classically, fMRI designs are optimized to maximize the direct effect. For example, the classic "block design" of alternating between two stimulus conditions every 30 seconds positions the expected neural modulation at an ideal temporal frequency with respect to the BOLD signal and noise properties.

To guide de Bruijn sequences for the detection of direct effects, construct the neural model to contain columns of expected relative magnitude of neural response. For example, an ever larger amplitude of response to the A, B, C, and D stimuli may be represented as:

| Syntax      | Description |
| ----------- | ----------- |
| Header      | Title       |
| Paragraph   | Text        |

^ ^A^B^C^D^
^A|1|2|3|4|
^B|1|2|3|4|
^C|1|2|3|4|
^D|1|2|3|4|

The resulting sequence will order the stimuli (within the constraints of counter-balance) to modulate the direct effect.

### Carry-over Effects
The transitions between stimuli may be expected to modulate neural response. These are *carry-over effects//((GK Aguirre (2007) [[http://www.ncbi.nlm.nih.gov/pubmed/17376705|Continuous carry-over designs for fMRI.]] *Neuroimage*, 35, 1480-1494.)). Neural adaptation (or habituation) is an example of a (symmetric) carry-over effect.

To design a sequence optimized for these effects, create a neural model which represents the effect of pair-wise transitions between the stimuli. The row index indicates the prior stimulus, and the column index indicates the current stimulus. This model, for example, contains both a symmetric (adaptation) effect, and an asymmetric bias effect:

^ ^A^B^C^D^
^A|0|1|2|3|
^B|1.5|0|1|2|
^C|2.5|1.5|0|1|
^D|3.5|2.5|1.5|0|

In this example, the carry over effects are predicted to be larger for the transition of ''[A ⇒ B]'' as compared to ''[B ⇒ A]'' (for example, transitions towards the "A" end of a linear stimulus space are more salient and perceived as larger).  The transition values represented are:

^prior stimulus^current stimulus^modeled transition^
|A|C|2|
|C|A|2.5|

Transitions between identical stimuli may be modeled as an undefined transition (''-1'') or as having a transition of zero (''0''), depending upon the hypotheses of the study and the modeling approach of the analysis.

### Guide function
Any valid de Bruijn cycle will provide the specified level of counter-balance amongst the labels, and thus stimuli, of the experiment. Not all orderings of stimuli, however, are equally useful for neuroimaging experiments. Because of the signal and noise properties of fMRI, some temporal frequencies of neural modulation are preferentially detected by the method. The *path-guided* approach to de Bruijn cycle generation encodes in the ordering of the stimuli a hypothesized neural modulation at these preferred frequencies.

To do so, we define a *guide function* in one of three ways:

  * Enter the word ''HRF''. A guide function will be internally generated that has the same power spectrum as the BOLD hemodynamic response function, although with no power below 0.01 Hz. This will generally produce a sequence with good detection power and a stochastic variation in stimulus transitions, although further improvements can be obtained by using a guide function that is positioned solely at a few or a single frequency. **NOTE** An SOA must be specified using the ''-eval'' option for this guide function to be defined.\\  \\  **-OR-**\\  \\
  * Enter a range of periods (in units of labels), as described by [Tmin,Tmax]. A guide function will be internally generated as the sum of sinusoids of random phase, each having a period equal to an integer between T<sub>min</sub> and T<sub>max</sub>. Enter the same value for T<sub>min</sub> and T<sub>max</sub> to guide the modulation at a single frequency. Generally, when ''(1e5 / SOA in ms) > T > (1e4 / SOA in ms)'' encoded modulations will be within a detectable range for the BOLD system.\\  \\  **-OR-**\\  \\  
  * Provide the path to an external file, which contains k<sup>n</sup> space-separated floating point values. Each element will correspond to a *relative* transition in the output sequence. The elements in the guide function may be either positive or negative; the vector will be normalized.


For fMRI, the optimal range of temporal frequencies is ~0.01-0.1 Hz((E Zarahn, GK Aguirre, M D'Esposito. (1997). [[http://www.ncbi.nlm.nih.gov/sites/entrez?Db=pubmed&Cmd=ShowDetailView&TermToSearch=9345548&|Empirical analyses of BOLD fMRI statistics. I. Spatially unsmoothed data collected under null-hypothesis conditions.]] *Neuroimage*, 5, 179-197.)). Higher frequencies are attenuated by the dispersed hemodynamic response, while lower frequencies are lost within the pink (1///f//) noise of the system. Perhaps surprisingly, it is not the case that the best performing sequence is always given by the lowest temporal frequency, even ignoring the presence of 1///f* noise. For some neural models, the available sets of transitions more readily fit higher temporal frequencies. This seems to be true in particular for one-dimensional stimulus spaces with a small number of stimuli. We recommend a search over a range of guide functions using a [[:public:de_bruijn_software#shell_script|shell script]].


The interpretation of the path-guided sequence in Hz requires a specification of the stimulus-onset asynchrony, described next.

### Evaluation mode and SOA
Given a particular hemodynamic response function, and a model of the 1///f* noise of the BOLD fMRI system, *detection power//((Thomas Liu. (2004) [[http://www.ncbi.nlm.nih.gov/pubmed/14741677|Efficiency, power, and entropy in event-related fMRI with multiple trial types. Part II: design of experiments]]. *Neuroimage* 21: 401-13.)) is the proportion of neural variance which appears in the imaging signal. It is calculated as the variance of the hypothesized neural modulation proportional to sequential stimulus distance as the denominator, and the numerator as the variance of that modulation after passing through the BOLD fMRI system. In our implementation, the BOLD system is modeled using a standard, population averaged hemodynamic response((GK Aguirre, E Zarahn, M D'Esposito. (1998). [[http://www.ncbi.nlm.nih.gov/pubmed/9811554|The variability of human BOLD hemodynamic responses.]] *NeuroImage*, 8, 360-369.)) and the elevated 1///f* noise range by a 0.01 Hz, high-pass notch filter.

When the ''-eval'' flag is set, the stimulus-onset asynchrony (SOA) parameter is specified in units of milliseconds (i.e., the time that elapses between the start of one stimulus and the start of the next stimulus). Along with the path-guided de Bruijn sequence, the routine then returns both the detection power, and the correlation coefficient of the guide function (input) with the sequence of distances between stimuli generated (output).

# Basic Examples

^Command^Interpretation^
|<code>./debruijn -h</code>|Displays the usage information.|
|<code>./debruijn 17 3</code>|Generates a deBruijn sequence with 17 labels and 3rd-level counterbalancing.|
|<code>./debruijn -v 17 3</code>|Generates a deBruijn sequence with 17 labels and 3rd-level counterbalancing, and prints output in verbose mode.|
|<code>./debruijn -t 10 2 5  [34,55] my_neural_model_matrix.txt</code>|Generates a deBruijn sequence with 10 labels and 2nd-level counterbalancing, and prints output in terse mode. The sequence is generated using 5 bins, a neural model matrix specified in the file *my_neural_model_matrix.txt*, and a guide function that is a sum of sinusoids with periods varying from 34 to 55 elements.|
|<code>./debruijn -t 10 2 5 HRF my_neural_model_matrix.txt -eval 1500</code>|Generates a deBruijn sequence with 10 labels and 2nd-level counterbalancing, and prints output in terse mode. The sequence is generated using 5 bins, a neural model matrix specified in the file *my_neural_model_matrix.txt*, and a guide function informed by the filtering properties of the BOLD hemodynamic response function is used. A stimulus-onset asynchrony of 1000 milliseconds is used in the evaluation of the sequences, and the theoretical detection power is returned.|
|<code>./debruijn 17 2 10 my_guide_function.txt my_neural_model_matrix.txt -eval 1000</code>|Generates a deBruijn sequence with 17 labels and 2nd-level counterbalancing, and prints output in normal mode. The sequence is generated using 10 bins, a neural model matrix specified in the file *my_neural_model_matrix.txt*, and a guide function specified in the file *my_guide_function.txt*. A stimulus-onset asynchrony of 1000 milliseconds is used in the evaluation of the sequences, and the theoretical detection power is returned.|

  ====== Example neural model matrices ======

===== Carry-over experiments =====
^  ^k=17, 16 stimuli di-oct^k=6, 5 stimuli linear^k=9, 8 element circular^
^Description|16 stimuli in a di-octagon arrangement (Euclidean geometry), stimulus zero as a null-trial, perfect repetitions have undefined distance|5 stimuli in a linear array, stimulus zero as null-trial, perfect repetitions have zero distance|8 stimuli with a circular similarity, stimulus zero as a null-trial, perfect repetitions have zero distance|
^Link to file|[[https://cfn.upenn.edu/aguirre/public/debruijn/k=17_16stim_dioct.txt|k=17_16stim_dioct.txt]]|[[https://cfn.upenn.edu/aguirre/public/debruijn/k=6_5stim_linear.txt|k=6_5stim_linear.txt]]|[[https://cfn.upenn.edu/aguirre/public/debruijn/k=9_8stim_circular.txt|k=9_8stim_circular.txt]]|
^Matrix|  {{:public:seqs:dioctsim.png?100}}  |  {{:public:seqs:5x5_lin_simspace.png?100}}  |  {{:public:seqs:8x8_circ_simspace.png?100}}  |
^Example|  {{:public:seqs:dioctspace.png?100}}  |  {{:public:seqs:5facelinearmorph.png?200}}  |  {{:public:seqs:8circ.png?200}}  |
^Reference|DM Drucker, WT Kerr, GK Aguirre. (2009) [[http://www.ncbi.nlm.nih.gov/pubmed/19357342?dopt=Abstract |Distinguishing conjoint and independent neural tuning for stimulus features with fMRI adaptation.]] *Journal of Neurophysiology*. June;101(6):3310-24|DA Kahn, AM Harris, DA Wolk, GK Aguirre. (2010) [[http://www.journalofvision.org/content/10/10/12.abstract |Temporally distinct neural coding of perceptual similarity and prototype bias.]] *Journal of Vision*, 10(10):12, 1-12.|  |


====== Shell script ======
This shell script may be used to search for and retain the best debruijn sequence, or to explore the effect of changes over the parameters of the search.

<code=bash>
#######################################
# Finds the best de Bruijn sequence
#######################################

clear

#######################################
# PARAMETERS:
#######################################
k=11
n=3
numbins=11
distfile=k=11_10stim_linear.txt
isi=1100
tmax=80
numiter=100
#######################################

maxdetecpow=0

for (( tmin = 10; tmin <= 80; tmin++))
do

for (( i = 1; i <= $numiter; i++ ))
do
	echo i=$i
        ./debruijn -t $k $n $numbins $distfile [$tmin,$tmax] -eval $isi -debug > myoutput.txt

	detecpow=`cat myoutput.txt | grep DETECTION | awk '{ print $3 }'`        #saves the detection power from the output text in a variable
	correlation=`cat myoutput.txt | grep CORRELATION | awk '{ print $3}'`    #saves the correlation from the output text in a variable

        if [ $detecpow ]                                                         #if the detection power was succesfully saved in a variable
        then
	        detecpow=`echo $detecpow | bc`                                   #converts from str to float
                compare_result=`echo "$detecpow > $maxdetecpow" | bc`            #compares with the current max value
                if test $compare_result -gt 0                                    #if the current one is greater
	        then
		        correlation=`echo $correlation | bc`                     #converts from str to float
                        echo CORRELATION = $correlation                          #outputs a text showing the correlation value for the sequence with max detection power
			echo DETECTION POWER = $detecpow                         #outputs a text showing the maximum detection power so far
			
			maxdetecpow=$detecpow                                    #saves the maximum detection power in a variable for future comparisons
			cp myoutput.txt bestoutput.txt                           #saves the output with the best detection power in bestoutput.txt
                fi
        fi

done

done
</code>
