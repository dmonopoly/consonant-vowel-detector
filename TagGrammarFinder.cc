// Main method outputs results (quick test).

#include "TagGrammarFinder.h"

#define DO_MAIN false
#define EXTRA_PRINTING false

bool TagGrammarFinder::GetTagGrammarFromOrganizedRows(
    const string &filename, map<Notation, double> *data,
    vector<string> *tag_list) {
  ifstream fin(filename.c_str());
  if (fin.fail()) {
    cerr << "Could not open file " << filename << endl;
    return false;
  } else {
    set<string> sounds;
    map<Notation, int> unigram_counts;
    map<Notation, int> bigram_counts;
    // Read bigram counts from file.
    int count;
    string sound1, sound2;
    // SUM_i C(s_i)
    Notation n_count_total(SIGMA + "C", {ARB_SOUND_PLACEHOLDER});
    while (true) {
      fin >> count;
      if (fin.eof())
        break;
      fin >> sound1 >> sound2;
      if (EXTRA_PRINTING)
        cout << "Read " << count << " " << sound1 << " " << sound2 << endl;
      sounds.insert(sound1);
      sounds.insert(sound2);
      Notation n_count_seq("C", {sound1}, TagGrammarFinder::SEQ_DELIM,
                           {sound2});
      Notation n_count_seq_total(SIGMA + "C", {sound1},
                                 TagGrammarFinder::SEQ_DELIM,
                                 {ARB_SOUND_PLACEHOLDER});
      bigram_counts[n_count_seq] = count; // should only encounter once 
      bigram_counts[n_count_seq_total] += count;
      // Single probabilities. Treat C(s1) = SUM_i C(s1 s_i).
      Notation n_count("C", {sound1});  // C(s1), followed by SUM_i C(s_i).
      unigram_counts[n_count] += count;
      unigram_counts[n_count_total] += count;
    }
    fin.close();
    // Determine tag grammar probabilities.
    // Unigrams.
    for (auto s1 = sounds.begin(); s1 != sounds.end(); ++s1) {
      Notation nSingle("P", {*s1});
      Notation n_count("C", {*s1});
      try {
        (*data)[nSingle] = (double) unigram_counts.at(n_count) / 
                                  unigram_counts.at(n_count_total);
      } catch (out_of_range &e) {
        cerr << "Out of range error in unigram stuff: " << e.what() << endl;
      }
    }
    // Bigrams.
    for (auto s1 = sounds.begin(); s1 != sounds.end(); ++s1) {
      for (auto s2 = sounds.begin(); s2 != sounds.end(); ++s2) {
        Notation nGiven("P", {*s2}, TagGrammarFinder::GIVEN_DELIM, {*s1});
        Notation n_count_seq("C", {*s1}, TagGrammarFinder::SEQ_DELIM, {*s2});
        Notation n_count_seq_total(SIGMA + "C", {*s1},
            TagGrammarFinder::SEQ_DELIM, {ARB_SOUND_PLACEHOLDER});
        // If no key found, then just set to 0.
        if (bigram_counts.find(n_count_seq) == bigram_counts.end() ||
            bigram_counts.find(n_count_seq_total) == bigram_counts.end()) {
          if (EXTRA_PRINTING) {
            if (bigram_counts.find(n_count_seq) == bigram_counts.end())
              cout << "Not found: " << n_count_seq << endl;
            if (bigram_counts.find(n_count_seq_total) ==
                bigram_counts.end())
              cout << "Not found: " << n_count_seq_total << endl;
          }
          (*data)[nGiven] = 0;
        } else {
          if (EXTRA_PRINTING) {
            cout << bigram_counts.at(n_count_seq) << "/" <<
              bigram_counts.at(n_count_seq_total) << endl;
          }
          double val = (double) bigram_counts.at(n_count_seq) /
            bigram_counts.at(n_count_seq_total);
          (*data)[nGiven] = val;
        }
      }
    }
    // Pass sounds to tag_list.
    for (auto it = sounds.begin(); it != sounds.end(); ++it) {
      tag_list->push_back(*it);
    }
  }
  return true;
}

// Quick test.
#if DO_MAIN
int main(int argc, char *argv[]) {
  if (argc < 2) {
    cerr << "No filename given." << endl;
    return 0;
  } 
  string filename = argv[1];
  map<Notation, double> data;
  vector<string> tag_list;
  TagGrammarFinder::GetTagGrammarFromOrganizedRows(filename, &data, &tag_list);
  cout << "Data:\n";
  for (auto it = data.begin(); it != data.end(); ++it) {
    cout << it->first << " " << it->second << endl;
  }
  cout << "Tag list:\n";
  for (auto it = tag_list.begin(); it != tag_list.end(); ++it) {
    cout << *it << ",";
  }
  cout << endl;
  return 0;
}
#endif
