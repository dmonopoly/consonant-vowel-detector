## What is this?
A program that distinguishes consonants vs. vowels in text.

### Purpose
Research. 

### Code Organization
- EMViterbiPackage contains the main code for running EM and Viterbi.
- Cyphers contains the text files that have the characters that we want to match
  Consonant and Vowel tags to.

## Running
Compile using cmake. Example:

    mkdir build && cd build
    cmake ..
    make

Then run on one of the cyphers in Cyphers/:

    ./consonant-vowel-detector ../Cyphers/eng.cypher.very_short.txt

### Options in Main.cc
At the top of Main.cc you can change the #defines to suit your purposes.

The following settings are sufficient to get nice results if you run on the
input file *corpus.spanish.quixote.written.txt*. This run, on my machine, took
about 1hr 20min, but the yielded results are a great matching of vowels and
consonants.

    #define EXTRA_PRINTING false  // Print extra things to debug.
    #define SHOW_PROBS_BEFORE_EM false
      // Print the probabilities set initially, prior to the first run of EM.
    #define WRITE_VITERBI_RESULTS_TO_FILE true
      // Write results to observed_data_probabilities.txt.
      // These probabilities should increase with each iteration.
    #define WRITE_LEARNED_PROBABILITIES true
      // Write final learned probabilities of the best run to
      // learned_probabilities_for_best_run.txt
    
    #define TREAT_UNDERSCORE_AS_SPACE true
      // If true, underscores in the cypher are marked as space tags (_') 100%
      // of the time. If false, this will try to identify what spaces are in the
      // cypher.
    #define NUMBER_ITERATIONS 20
    #define RANDOM_INITIAL_START true  // false means uniform probs used.
    #define NUM_RESTARTS 20  // Used only if RANDOM_INITIAL_START is true.
    #define PRINT_RESULTS_OF_EACH_RANDOM_RESTART true
      // Prints results to the terminal as the program runs. Handy so you can
      // see the progress of the program.

### Options in TrellisAid.cc
You can similarly change options in TrellisAid.cc, which has the actual
EM and Viterbi implementation. The following are typically fine though.

    #define EXTRA_PRINTING false
    #define PRINT_VITERBI_RESULTS_OFTEN false

You might change the latter to true if you want to see the specific results of
each viterbi run, like the result after 1 iteration of EM, then 2, etc. But
false is good if you have lots of random restarts and/or lots of iterations.

## Notes on results
The following tags are used in the system:

    C', V', and _'

- Underscore means space.
- The apostrophe is there to differentiate it from an observed symbol from the
  cypher. e.g., C is the observed symbol found from a cypher, while C with an
  apostrophe is the tag for Consonant.
- *The tags are arbitrary symbols,* so C' and V' may be mixed up - C' may match
  vowels and V' may match consonants. The program just finds the most probable
  matching.

Log probabilities are used, so expect those in the output files

    observed_data_probabilities.txt
    learned_probabilities_for_best_run.txt

From the latter file, we can make statements like "20% of V' are realized with
the letter 'a'" or "yes, everything adds to one" or "SPACE is followed by C 72%
of the time." Converting to regular probabilities (which can be done only if
probabilities aren't too small) would make seeing these statements easier.
