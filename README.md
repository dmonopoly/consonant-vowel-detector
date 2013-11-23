## What is this?
A program that distinguishes consonants vs. vowels in text.

### Purpose
Research. 

TODO: update all below.

### Code Organization
- EMViterbiPackage contains the main code for running EM and Viterbi.
- GeneratedData contains example outputs.
- LanguageModel contains count data.
- Cyphers contains cyphers to break.

## Running
Compile using cmake. Example:

    mkdir build && cd build
    cmake ..
    make

Run the generated executable. Alternatively run any test executables, but make
sure you do `cmake -Dtest=ON ..` first.

Running Main:

  Example 1 - Spanish cypher (all one line):

    ./sound-decipher ../LanguageModels/corpus.spanish.sound_bigrams.txt
    ../Cyphers/corpus.spanish.quixote.written_very_short.txt > output.txt

  Running on the very short sequence takes roughly 20 seconds. For long
  sequences, it can take quite a while.

  Example 2 - English cypher (all one line):

    ./sound-decipher ../LanguageModels/eng.bigram.counts.txt
    ../Cyphers/eng.cypher.short.txt > eng-output.txt

CypherReader can be compiled and run independently:

    # Uncomment main() function in CypherReader.cc. Then:
    clang++ -std=c++11 -stdlib=libc++ CypherReader.cc

TagGrammarFinder also:

    # First change "DO_MAIN false" to true in TagGrammarFinder.cc. Then:
    cmake .. -Dtag-grammar=on
    make
    ./tag-grammar-finder ../LanguageModels/eng.bigram.counts.txt

    # Revert back to normal by changing DO_MAIN back to true, and then
    cmake .. -Dtag-grammar=off
