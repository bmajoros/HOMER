CC 		= g++
LDFLAGS	        = -O
CFLAGS		= -O
OBJ		= obj
TIGRPLUSPLUS	= tigr++
HOMER		= .
LIBS		= 

all: \
	obj \
	deterministic-trainer \
	homer \
	baum-welch \
	model-combiner \
	perl-files \
	warning \

.PHONY : perl-files
perl-files:
	chmod +x *.pl

.PHONY : warning
warning:
	@echo WARNING!  Do not attempt to train the gene-finder in the installation directory.  Read the included README file.  An auto-trainer script is provided for your convenience.

obj:
	mkdir obj

.PHONY : clean
clean:
	@rm obj/*.o


$(OBJ)/TigrFastaReader.o: \
		$(TIGRPLUSPLUS)/TigrFastaReader.H \
		$(TIGRPLUSPLUS)/TigrFastaReader.C
	g++ $(CFLAGS) -o $(OBJ)/TigrFastaReader.o -c \
		$(TIGRPLUSPLUS)/TigrFastaReader.C

$(OBJ)/DnaAlphabet.o: \
		$(TIGRPLUSPLUS)/DnaAlphabet.H \
		$(TIGRPLUSPLUS)/DnaAlphabet.C
	g++ $(CFLAGS) -o $(OBJ)/DnaAlphabet.o -c \
		$(TIGRPLUSPLUS)/DnaAlphabet.C

$(OBJ)/NthOrderStringIterator.o: \
		$(TIGRPLUSPLUS)/NthOrderStringIterator.H \
		$(TIGRPLUSPLUS)/NthOrderStringIterator.C
	g++ $(CFLAGS) -o $(OBJ)/NthOrderStringIterator.o -c \
		$(TIGRPLUSPLUS)/NthOrderStringIterator.C

$(OBJ)/model-combiner.o: \
		$(HOMER)/model-combiner.C \
		$(TIGRPLUSPLUS)/Alphabet.H \
		$(HOMER)/ForwardAlgorithm.H \
		$(HOMER)/BaumWelch.H \
		$(HOMER)/HiddenMarkovModel.H \
		$(TIGRPLUSPLUS)/Sequence.H
	g++ $(CFLAGS) -o $(OBJ)/model-combiner.o -c \
		$(HOMER)/model-combiner.C

$(OBJ)/deterministic-trainer.o: \
		$(HOMER)/deterministic-trainer.C \
		$(TIGRPLUSPLUS)/Alphabet.H \
		$(HOMER)/ForwardAlgorithm.H \
		$(HOMER)/BaumWelch.H \
		$(HOMER)/HiddenMarkovModel.H \
		$(TIGRPLUSPLUS)/Sequence.H
	g++ $(CFLAGS) -o $(OBJ)/deterministic-trainer.o -c \
		$(HOMER)/deterministic-trainer.C

$(OBJ)/baum-welch.o: \
		$(HOMER)/baum-welch.C \
		$(TIGRPLUSPLUS)/Alphabet.H \
		$(HOMER)/ForwardAlgorithm.H \
		$(HOMER)/BaumWelch.H \
		$(HOMER)/HiddenMarkovModel.H \
		$(TIGRPLUSPLUS)/Sequence.H
	g++ $(CFLAGS) -o $(OBJ)/baum-welch.o -c \
		$(HOMER)/baum-welch.C

$(OBJ)/MCE.o: \
		$(HOMER)/MCE.C \
		$(TIGRPLUSPLUS)/Alphabet.H \
		$(HOMER)/HiddenMarkovModel.H \
		$(TIGRPLUSPLUS)/Sequence.H
	g++ $(CFLAGS) -o $(OBJ)/MCE.o -c \
		$(HOMER)/MCE.C

$(OBJ)/hmm-substitute.o: \
		$(HOMER)/hmm-substitute.C \
		$(TIGRPLUSPLUS)/Alphabet.H \
		$(HOMER)/HiddenMarkovModel.H \
		$(TIGRPLUSPLUS)/Sequence.H
	g++ $(CFLAGS) -o $(OBJ)/hmm-substitute.o -c \
		$(HOMER)/hmm-substitute.C

$(OBJ)/Alphabet.o: \
		$(TIGRPLUSPLUS)/Alphabet.H \
		$(TIGRPLUSPLUS)/Alphabet.C
	g++ $(CFLAGS) -o $(OBJ)/Alphabet.o -c \
		$(TIGRPLUSPLUS)/Alphabet.C

$(OBJ)/HMMbuilder.o: \
		$(HOMER)/HMMbuilder.H \
		$(HOMER)/HMMbuilder.C
	g++ $(CFLAGS) -o $(OBJ)/HMMbuilder.o -c \
		$(HOMER)/HMMbuilder.C

$(OBJ)/HigherOrderAlphabet.o: \
		$(HOMER)/HigherOrderAlphabet.H \
		$(HOMER)/HigherOrderAlphabet.C
	g++ $(CFLAGS) -o $(OBJ)/HigherOrderAlphabet.o -c \
		$(HOMER)/HigherOrderAlphabet.C

$(OBJ)/HMMreader.o: \
		$(HOMER)/HMMreader.H \
		$(HOMER)/HMMreader.C
	g++ $(CFLAGS) -o $(OBJ)/HMMreader.o -c \
		$(HOMER)/HMMreader.C

$(OBJ)/BackwardAlgorithm.o: \
		$(HOMER)/HiddenMarkovModel.H \
		$(TIGRPLUSPLUS)/Sequence.H \
		$(HOMER)/BackwardAlgorithm.H \
		$(HOMER)/BackwardAlgorithm.C
	g++ $(CFLAGS) -o $(OBJ)/BackwardAlgorithm.o -c \
		$(HOMER)/BackwardAlgorithm.C

$(OBJ)/BaumWelch.o: \
		$(HOMER)/HiddenMarkovModel.H \
		$(TIGRPLUSPLUS)/Sequence.H \
		$(HOMER)/ForwardAlgorithm.H \
		$(HOMER)/BackwardAlgorithm.H \
		$(HOMER)/BaumWelch.H \
		$(HOMER)/BaumWelch.C
	g++ $(CFLAGS) -o $(OBJ)/BaumWelch.o -c \
		$(HOMER)/BaumWelch.C

$(OBJ)/ForwardAlgorithm.o: \
		$(TIGRPLUSPLUS)/Sequence.H \
		$(HOMER)/HiddenMarkovModel.H \
		$(HOMER)/ForwardAlgorithm.H \
		$(HOMER)/ForwardAlgorithm.C
	g++ $(CFLAGS) -o $(OBJ)/ForwardAlgorithm.o -c \
		$(HOMER)/ForwardAlgorithm.C

$(OBJ)/FastViterbi.o: \
		$(TIGRPLUSPLUS)/Sequence.H \
		$(HOMER)/HiddenMarkovModel.H \
		$(HOMER)/HMMGraph.H \
		$(HOMER)/FastViterbi.H \
		$(HOMER)/FastViterbi.C
	g++ $(CFLAGS) -o $(OBJ)/FastViterbi.o -c \
		$(HOMER)/FastViterbi.C

$(OBJ)/HiddenMarkovModel.o: \
		$(HOMER)/HiddenMarkovModel.H \
		$(HOMER)/HiddenMarkovModel.C \
		$(TIGRPLUSPLUS)/Alphabet.H \
		$(TIGRPLUSPLUS)/Symbol.H \
		$(TIGRPLUSPLUS)/TigrRandom.H
	g++ $(CFLAGS) -o $(OBJ)/HiddenMarkovModel.o -c \
		$(HOMER)/HiddenMarkovModel.C

$(OBJ)/Sequence.o: \
		$(TIGRPLUSPLUS)/Sequence.H \
		$(TIGRPLUSPLUS)/Symbol.H \
		$(TIGRPLUSPLUS)/Alphabet.H \
		$(TIGRPLUSPLUS)/Sequence.C
	g++ $(CFLAGS) -o $(OBJ)/Sequence.o -c \
		$(TIGRPLUSPLUS)/Sequence.C

$(OBJ)/Symbol.o: \
		$(TIGRPLUSPLUS)/Symbol.H \
		$(TIGRPLUSPLUS)/Symbol.C
	g++ $(CFLAGS) -o $(OBJ)/Symbol.o -c \
		$(TIGRPLUSPLUS)/Symbol.C

$(OBJ)/TigrRandom.o: \
		$(TIGRPLUSPLUS)/TigrRandom.H \
		$(TIGRPLUSPLUS)/TigrRandom.C
	g++ $(CFLAGS) -o $(OBJ)/TigrRandom.o -c \
		$(TIGRPLUSPLUS)/TigrRandom.C

$(OBJ)/TigrProgress.o: \
		$(TIGRPLUSPLUS)/TigrProgress.H \
		$(TIGRPLUSPLUS)/TigrProgress.C
	g++ $(CFLAGS) -o $(OBJ)/TigrProgress.o -c \
		$(TIGRPLUSPLUS)/TigrProgress.C

$(OBJ)/TigrCommandLine.o: \
		$(TIGRPLUSPLUS)/TigrCommandLine.H \
		$(TIGRPLUSPLUS)/TigrCommandLine.C
	g++ $(CFLAGS) -o $(OBJ)/TigrCommandLine.o -c \
		$(TIGRPLUSPLUS)/TigrCommandLine.C

$(OBJ)/TigrConfigFile.o: \
		$(TIGRPLUSPLUS)/TigrConfigFile.H \
		$(TIGRPLUSPLUS)/TigrConfigFile.C
	g++ $(CFLAGS) -o $(OBJ)/TigrConfigFile.o -c \
		$(TIGRPLUSPLUS)/TigrConfigFile.C

$(OBJ)/TigrStrTokenizer.o: \
		$(TIGRPLUSPLUS)/TigrStrTokenizer.H \
		$(TIGRPLUSPLUS)/TigrStrTokenizer.C
	g++ $(CFLAGS) -o $(OBJ)/TigrStrTokenizer.o -c \
		$(TIGRPLUSPLUS)/TigrStrTokenizer.C

$(OBJ)/TigrFile.o: \
		$(TIGRPLUSPLUS)/TigrFile.H \
		$(TIGRPLUSPLUS)/TigrFile.C
	g++ $(CFLAGS) -o $(OBJ)/TigrFile.o -c \
		$(TIGRPLUSPLUS)/TigrFile.C

$(OBJ)/TigrString.o: \
		$(TIGRPLUSPLUS)/TigrString.H \
		$(TIGRPLUSPLUS)/TigrString.C
	g++ $(CFLAGS) -o $(OBJ)/TigrString.o -c \
		$(TIGRPLUSPLUS)/TigrString.C

$(OBJ)/TigrRegex.o: \
		$(TIGRPLUSPLUS)/TigrRegex.H \
		$(TIGRPLUSPLUS)/TigrRegex.C
	g++ $(CFLAGS) -o $(OBJ)/TigrRegex.o -c \
		$(TIGRPLUSPLUS)/TigrRegex.C

$(OBJ)/TigrExceptions.o: \
		$(TIGRPLUSPLUS)/TigrExceptions.H \
		$(TIGRPLUSPLUS)/TigrExceptions.C
	g++ $(CFLAGS) -o $(OBJ)/TigrExceptions.o -c \
		$(TIGRPLUSPLUS)/TigrExceptions.C

$(OBJ)/TigrStacktrace.o: \
		$(TIGRPLUSPLUS)/TigrStacktrace.H \
		$(TIGRPLUSPLUS)/TigrStacktrace.C
	g++ $(CFLAGS) -o $(OBJ)/TigrStacktrace.o -c \
		$(TIGRPLUSPLUS)/TigrStacktrace.C

$(OBJ)/TigrBitSet.o: \
		$(TIGRPLUSPLUS)/TigrBitSet.H \
		$(TIGRPLUSPLUS)/TigrBitSet.C
	g++ $(CFLAGS) -o $(OBJ)/TigrBitSet.o -c \
		$(TIGRPLUSPLUS)/TigrBitSet.C

$(OBJ)/TigrGffReader.o: \
		$(TIGRPLUSPLUS)/TigrGffReader.H \
		$(TIGRPLUSPLUS)/TigrGffReader.C
	g++ $(CFLAGS) -o $(OBJ)/TigrGffReader.o -c \
		$(TIGRPLUSPLUS)/TigrGffReader.C
$(OBJ)/TigrGffExon.o: \
		$(TIGRPLUSPLUS)/TigrGffExon.H \
		$(TIGRPLUSPLUS)/TigrGffExon.C
	g++ $(CFLAGS) -o $(OBJ)/TigrGffExon.o -c \
		$(TIGRPLUSPLUS)/TigrGffExon.C
$(OBJ)/TigrGffTranscript.o: \
		$(TIGRPLUSPLUS)/TigrGffTranscript.H \
		$(TIGRPLUSPLUS)/TigrGffTranscript.C
	g++ $(CFLAGS) -o $(OBJ)/TigrGffTranscript.o -c \
		$(TIGRPLUSPLUS)/TigrGffTranscript.C
$(OBJ)/TigrGffFeature.o: \
		$(TIGRPLUSPLUS)/TigrGffFeature.H \
		$(TIGRPLUSPLUS)/TigrGffFeature.C
	g++ $(CFLAGS) -o $(OBJ)/TigrGffFeature.o -c \
		$(TIGRPLUSPLUS)/TigrGffFeature.C

baum-welch: \
		$(OBJ)/DnaAlphabet.o \
		$(OBJ)/NthOrderStringIterator.o \
		$(OBJ)/HigherOrderAlphabet.o \
		$(OBJ)/HMMGraph.o \
		$(OBJ)/FastForward.o \
		$(OBJ)/baum-welch.o \
		$(OBJ)/TigrFastaReader.o \
		$(OBJ)/TigrRegex.o \
		$(OBJ)/Alphabet.o \
		$(OBJ)/BackwardAlgorithm.o \
		$(OBJ)/BaumWelch.o \
		$(OBJ)/ForwardAlgorithm.o \
		$(OBJ)/HiddenMarkovModel.o \
		$(OBJ)/Sequence.o \
		$(OBJ)/TigrProgress.o \
		$(OBJ)/TigrRandom.o \
		$(OBJ)/TigrStrTokenizer.o \
		$(OBJ)/TigrFile.o \
		$(OBJ)/TigrConfigFile.o \
		$(OBJ)/TigrString.o \
		$(OBJ)/TigrStacktrace.o \
		$(OBJ)/TigrCommandLine.o \
		$(OBJ)/HMMbuilder.o \
		$(OBJ)/HMMreader.o \
		$(OBJ)/Symbol.o
	g++ $(LDFLAGS) -o baum-welch \
		$(OBJ)/HMMGraph.o \
		$(OBJ)/NthOrderStringIterator.o \
		$(OBJ)/HigherOrderAlphabet.o \
		$(OBJ)/FastForward.o \
		$(OBJ)/baum-welch.o \
		$(OBJ)/TigrFastaReader.o \
		$(OBJ)/TigrRegex.o \
		$(OBJ)/TigrString.o \
		$(OBJ)/HMMreader.o \
		$(OBJ)/HMMbuilder.o \
		$(OBJ)/TigrCommandLine.o \
		$(OBJ)/TigrProgress.o \
		$(OBJ)/Alphabet.o \
		$(OBJ)/TigrConfigFile.o \
		$(OBJ)/BackwardAlgorithm.o \
		$(OBJ)/BaumWelch.o \
		$(OBJ)/ForwardAlgorithm.o \
		$(OBJ)/DnaAlphabet.o \
		$(OBJ)/TigrRandom.o \
		$(OBJ)/TigrStacktrace.o \
		$(OBJ)/TigrStrTokenizer.o \
		$(OBJ)/TigrFile.o \
		$(OBJ)/HiddenMarkovModel.o \
		$(OBJ)/Sequence.o \
		$(OBJ)/Symbol.o \
		$(LIBS)

MCE: \
		$(OBJ)/FastViterbi.o \
		$(OBJ)/DnaAlphabet.o \
		$(OBJ)/NthOrderStringIterator.o \
		$(OBJ)/HigherOrderAlphabet.o \
		$(OBJ)/HMMGraph.o \
		$(OBJ)/FastForward.o \
		$(OBJ)/MCE.o \
		$(OBJ)/TigrFastaReader.o \
		$(OBJ)/TigrRegex.o \
		$(OBJ)/Alphabet.o \
		$(OBJ)/BackwardAlgorithm.o \
		$(OBJ)/BaumWelch.o \
		$(OBJ)/ForwardAlgorithm.o \
		$(OBJ)/HiddenMarkovModel.o \
		$(OBJ)/Sequence.o \
		$(OBJ)/TigrProgress.o \
		$(OBJ)/TigrRandom.o \
		$(OBJ)/TigrStrTokenizer.o \
		$(OBJ)/TigrFile.o \
		$(OBJ)/TigrConfigFile.o \
		$(OBJ)/TigrString.o \
		$(OBJ)/TigrStacktrace.o \
		$(OBJ)/TigrCommandLine.o \
		$(OBJ)/HMMbuilder.o \
		$(OBJ)/HMMreader.o \
		$(OBJ)/Symbol.o
	g++ $(LDFLAGS) -o MCE \
		$(OBJ)/FastViterbi.o \
		$(OBJ)/HMMGraph.o \
		$(OBJ)/NthOrderStringIterator.o \
		$(OBJ)/HigherOrderAlphabet.o \
		$(OBJ)/FastForward.o \
		$(OBJ)/MCE.o \
		$(OBJ)/TigrFastaReader.o \
		$(OBJ)/TigrRegex.o \
		$(OBJ)/TigrString.o \
		$(OBJ)/HMMreader.o \
		$(OBJ)/HMMbuilder.o \
		$(OBJ)/TigrCommandLine.o \
		$(OBJ)/TigrProgress.o \
		$(OBJ)/Alphabet.o \
		$(OBJ)/TigrConfigFile.o \
		$(OBJ)/BackwardAlgorithm.o \
		$(OBJ)/BaumWelch.o \
		$(OBJ)/ForwardAlgorithm.o \
		$(OBJ)/DnaAlphabet.o \
		$(OBJ)/TigrRandom.o \
		$(OBJ)/TigrStacktrace.o \
		$(OBJ)/TigrStrTokenizer.o \
		$(OBJ)/TigrFile.o \
		$(OBJ)/HiddenMarkovModel.o \
		$(OBJ)/Sequence.o \
		$(OBJ)/Symbol.o \
		$(LIBS)

hmm-substitute: \
		$(OBJ)/FastViterbi.o \
		$(OBJ)/DnaAlphabet.o \
		$(OBJ)/NthOrderStringIterator.o \
		$(OBJ)/HigherOrderAlphabet.o \
		$(OBJ)/HMMGraph.o \
		$(OBJ)/FastForward.o \
		$(OBJ)/hmm-substitute.o \
		$(OBJ)/TigrFastaReader.o \
		$(OBJ)/TigrRegex.o \
		$(OBJ)/Alphabet.o \
		$(OBJ)/BackwardAlgorithm.o \
		$(OBJ)/BaumWelch.o \
		$(OBJ)/ForwardAlgorithm.o \
		$(OBJ)/HiddenMarkovModel.o \
		$(OBJ)/Sequence.o \
		$(OBJ)/TigrProgress.o \
		$(OBJ)/TigrRandom.o \
		$(OBJ)/TigrStrTokenizer.o \
		$(OBJ)/TigrFile.o \
		$(OBJ)/TigrConfigFile.o \
		$(OBJ)/TigrString.o \
		$(OBJ)/TigrStacktrace.o \
		$(OBJ)/TigrCommandLine.o \
		$(OBJ)/HMMbuilder.o \
		$(OBJ)/HMMreader.o \
		$(OBJ)/Symbol.o
	g++ $(LDFLAGS) -o hmm-substitute \
		$(OBJ)/FastViterbi.o \
		$(OBJ)/HMMGraph.o \
		$(OBJ)/NthOrderStringIterator.o \
		$(OBJ)/HigherOrderAlphabet.o \
		$(OBJ)/FastForward.o \
		$(OBJ)/hmm-substitute.o \
		$(OBJ)/TigrFastaReader.o \
		$(OBJ)/TigrRegex.o \
		$(OBJ)/TigrString.o \
		$(OBJ)/HMMreader.o \
		$(OBJ)/HMMbuilder.o \
		$(OBJ)/TigrCommandLine.o \
		$(OBJ)/TigrProgress.o \
		$(OBJ)/Alphabet.o \
		$(OBJ)/TigrConfigFile.o \
		$(OBJ)/BackwardAlgorithm.o \
		$(OBJ)/BaumWelch.o \
		$(OBJ)/ForwardAlgorithm.o \
		$(OBJ)/DnaAlphabet.o \
		$(OBJ)/TigrRandom.o \
		$(OBJ)/TigrStacktrace.o \
		$(OBJ)/TigrStrTokenizer.o \
		$(OBJ)/TigrFile.o \
		$(OBJ)/HiddenMarkovModel.o \
		$(OBJ)/Sequence.o \
		$(OBJ)/Symbol.o \
		$(LIBS)

deterministic-trainer: \
		$(OBJ)/HMMGraph.o \
		$(OBJ)/FastForward.o \
		$(OBJ)/deterministic-trainer.o \
		$(OBJ)/TigrFastaReader.o \
		$(OBJ)/TigrRegex.o \
		$(OBJ)/Alphabet.o \
		$(OBJ)/BackwardAlgorithm.o \
		$(OBJ)/BaumWelch.o \
		$(OBJ)/ForwardAlgorithm.o \
		$(OBJ)/DnaAlphabet.o \
		$(OBJ)/HiddenMarkovModel.o \
		$(OBJ)/Sequence.o \
		$(OBJ)/TigrProgress.o \
		$(OBJ)/TigrRandom.o \
		$(OBJ)/TigrStrTokenizer.o \
		$(OBJ)/TigrFile.o \
		$(OBJ)/TigrConfigFile.o \
		$(OBJ)/TigrString.o \
		$(OBJ)/TigrStacktrace.o \
		$(OBJ)/TigrCommandLine.o \
		$(OBJ)/HMMbuilder.o \
		$(OBJ)/HMMreader.o \
		$(OBJ)/Symbol.o
	g++ $(LDFLAGS) -o deterministic-trainer \
		$(OBJ)/HMMGraph.o \
		$(OBJ)/FastForward.o \
		$(OBJ)/deterministic-trainer.o \
		$(OBJ)/TigrFastaReader.o \
		$(OBJ)/TigrRegex.o \
		$(OBJ)/TigrString.o \
		$(OBJ)/HMMreader.o \
		$(OBJ)/HMMbuilder.o \
		$(OBJ)/TigrCommandLine.o \
		$(OBJ)/TigrProgress.o \
		$(OBJ)/Alphabet.o \
		$(OBJ)/TigrConfigFile.o \
		$(OBJ)/DnaAlphabet.o \
		$(OBJ)/BackwardAlgorithm.o \
		$(OBJ)/BaumWelch.o \
		$(OBJ)/ForwardAlgorithm.o \
		$(OBJ)/TigrRandom.o \
		$(OBJ)/TigrStacktrace.o \
		$(OBJ)/TigrStrTokenizer.o \
		$(OBJ)/TigrFile.o \
		$(OBJ)/HiddenMarkovModel.o \
		$(OBJ)/Sequence.o \
		$(OBJ)/Symbol.o \
		$(LIBS)

model-combiner: \
		$(OBJ)/HMMGraph.o \
		$(OBJ)/FastForward.o \
		$(OBJ)/model-combiner.o \
		$(OBJ)/TigrFastaReader.o \
		$(OBJ)/TigrRegex.o \
		$(OBJ)/Alphabet.o \
		$(OBJ)/BackwardAlgorithm.o \
		$(OBJ)/BaumWelch.o \
		$(OBJ)/ForwardAlgorithm.o \
		$(OBJ)/HiddenMarkovModel.o \
		$(OBJ)/Sequence.o \
		$(OBJ)/TigrProgress.o \
		$(OBJ)/TigrRandom.o \
		$(OBJ)/TigrStrTokenizer.o \
		$(OBJ)/TigrFile.o \
		$(OBJ)/DnaAlphabet.o \
		$(OBJ)/TigrConfigFile.o \
		$(OBJ)/TigrString.o \
		$(OBJ)/TigrStacktrace.o \
		$(OBJ)/TigrCommandLine.o \
		$(OBJ)/HMMbuilder.o \
		$(OBJ)/HMMreader.o \
		$(OBJ)/Symbol.o
	g++ $(LDFLAGS) -o model-combiner \
		$(OBJ)/HMMGraph.o \
		$(OBJ)/FastForward.o \
		$(OBJ)/model-combiner.o \
		$(OBJ)/TigrFastaReader.o \
		$(OBJ)/TigrRegex.o \
		$(OBJ)/TigrString.o \
		$(OBJ)/HMMreader.o \
		$(OBJ)/HMMbuilder.o \
		$(OBJ)/TigrCommandLine.o \
		$(OBJ)/TigrProgress.o \
		$(OBJ)/Alphabet.o \
		$(OBJ)/TigrConfigFile.o \
		$(OBJ)/BackwardAlgorithm.o \
		$(OBJ)/BaumWelch.o \
		$(OBJ)/DnaAlphabet.o \
		$(OBJ)/ForwardAlgorithm.o \
		$(OBJ)/TigrRandom.o \
		$(OBJ)/TigrStacktrace.o \
		$(OBJ)/TigrStrTokenizer.o \
		$(OBJ)/TigrFile.o \
		$(OBJ)/HiddenMarkovModel.o \
		$(OBJ)/Sequence.o \
		$(OBJ)/Symbol.o \
		$(LIBS)
#--------------------------------------------------------
$(OBJ)/homer.o:\
		$(HOMER)/homer.C
	$(CC) $(CFLAGS) -o $(OBJ)/homer.o -c \
		$(HOMER)/homer.C
#--------------------------------------------------------
$(OBJ)/score-best-path.o:\
		$(HOMER)/score-best-path.C
	$(CC) $(CFLAGS) -o $(OBJ)/score-best-path.o -c \
		$(HOMER)/score-best-path.C
#--------------------------------------------------------
$(OBJ)/HMMGraph.o:\
		$(HOMER)/HMMGraph.C \
		$(HOMER)/HMMGraph.H
	$(CC) $(CFLAGS) -o $(OBJ)/HMMGraph.o -c \
		$(HOMER)/HMMGraph.C
#--------------------------------------------------------
$(OBJ)/FastForward.o:\
		$(HOMER)/FastForward.C \
		$(HOMER)/FastForward.H
	$(CC) $(CFLAGS) -o $(OBJ)/FastForward.o -c \
		$(HOMER)/FastForward.C
#---------------------------------------------------------
homer: \
		$(OBJ)/HiddenMarkovModel.o \
		$(OBJ)/HigherOrderAlphabet.o \
		$(OBJ)/NthOrderStringIterator.o \
		$(OBJ)/HMMGraph.o \
		$(OBJ)/FastViterbi.o \
		$(OBJ)/FastForward.o \
		$(OBJ)/HMMreader.o \
		$(OBJ)/HMMbuilder.o \
		$(OBJ)/TigrGffReader.o \
		$(OBJ)/TigrGffExon.o \
		$(OBJ)/TigrGffFeature.o \
		$(OBJ)/TigrGffTranscript.o \
		$(OBJ)/TigrExceptions.o \
		$(OBJ)/TigrBitSet.o \
		$(OBJ)/TigrCommandLine.o \
		$(OBJ)/TigrString.o \
		$(OBJ)/homer.o \
		$(OBJ)/TigrFastaReader.o \
		$(OBJ)/DnaAlphabet.o \
		$(OBJ)/TigrRegex.o \
		$(OBJ)/Alphabet.o \
		$(OBJ)/Sequence.o \
		$(OBJ)/TigrProgress.o \
		$(OBJ)/TigrRandom.o \
		$(OBJ)/TigrStrTokenizer.o \
		$(OBJ)/TigrFile.o \
		$(OBJ)/TigrConfigFile.o \
		$(OBJ)/TigrStacktrace.o \
		$(OBJ)/Symbol.o
	$(CC) $(LDFLAGS) -o homer \
		$(OBJ)/TigrGffReader.o \
		$(OBJ)/HigherOrderAlphabet.o \
		$(OBJ)/NthOrderStringIterator.o \
		$(OBJ)/FastForward.o \
		$(OBJ)/TigrGffExon.o \
		$(OBJ)/TigrGffFeature.o \
		$(OBJ)/TigrGffTranscript.o \
		$(OBJ)/TigrExceptions.o \
		$(OBJ)/TigrBitSet.o \
		$(OBJ)/TigrCommandLine.o \
		$(OBJ)/TigrString.o \
		$(OBJ)/homer.o \
		$(OBJ)/TigrFastaReader.o \
		$(OBJ)/DnaAlphabet.o \
		$(OBJ)/TigrRegex.o \
		$(OBJ)/Alphabet.o \
		$(OBJ)/HiddenMarkovModel.o \
		$(OBJ)/HMMGraph.o \
		$(OBJ)/Sequence.o \
		$(OBJ)/TigrProgress.o \
		$(OBJ)/FastViterbi.o \
		$(OBJ)/TigrRandom.o \
		$(OBJ)/TigrStrTokenizer.o \
		$(OBJ)/TigrFile.o \
		$(OBJ)/TigrConfigFile.o \
		$(OBJ)/TigrStacktrace.o \
		$(OBJ)/HMMbuilder.o \
		$(OBJ)/HMMreader.o \
		$(OBJ)/Symbol.o \
		$(LIBS)
#--------------------------------------------------------
#---------------------------------------------------------
score-best-path: \
		$(OBJ)/HiddenMarkovModel.o \
		$(OBJ)/HigherOrderAlphabet.o \
		$(OBJ)/NthOrderStringIterator.o \
		$(OBJ)/HMMGraph.o \
		$(OBJ)/FastViterbi.o \
		$(OBJ)/FastForward.o \
		$(OBJ)/HMMreader.o \
		$(OBJ)/HMMbuilder.o \
		$(OBJ)/TigrGffReader.o \
		$(OBJ)/TigrGffExon.o \
		$(OBJ)/TigrGffFeature.o \
		$(OBJ)/TigrGffTranscript.o \
		$(OBJ)/TigrExceptions.o \
		$(OBJ)/TigrBitSet.o \
		$(OBJ)/TigrCommandLine.o \
		$(OBJ)/TigrString.o \
		$(OBJ)/score-best-path.o \
		$(OBJ)/TigrFastaReader.o \
		$(OBJ)/DnaAlphabet.o \
		$(OBJ)/TigrRegex.o \
		$(OBJ)/Alphabet.o \
		$(OBJ)/Sequence.o \
		$(OBJ)/TigrProgress.o \
		$(OBJ)/TigrRandom.o \
		$(OBJ)/TigrStrTokenizer.o \
		$(OBJ)/TigrFile.o \
		$(OBJ)/TigrConfigFile.o \
		$(OBJ)/TigrStacktrace.o \
		$(OBJ)/Symbol.o
	$(CC) $(LDFLAGS) -o score-best-path \
		$(OBJ)/TigrGffReader.o \
		$(OBJ)/HigherOrderAlphabet.o \
		$(OBJ)/NthOrderStringIterator.o \
		$(OBJ)/FastForward.o \
		$(OBJ)/TigrGffExon.o \
		$(OBJ)/TigrGffFeature.o \
		$(OBJ)/TigrGffTranscript.o \
		$(OBJ)/TigrExceptions.o \
		$(OBJ)/TigrBitSet.o \
		$(OBJ)/TigrCommandLine.o \
		$(OBJ)/TigrString.o \
		$(OBJ)/score-best-path.o \
		$(OBJ)/TigrFastaReader.o \
		$(OBJ)/DnaAlphabet.o \
		$(OBJ)/TigrRegex.o \
		$(OBJ)/Alphabet.o \
		$(OBJ)/HiddenMarkovModel.o \
		$(OBJ)/HMMGraph.o \
		$(OBJ)/Sequence.o \
		$(OBJ)/TigrProgress.o \
		$(OBJ)/FastViterbi.o \
		$(OBJ)/TigrRandom.o \
		$(OBJ)/TigrStrTokenizer.o \
		$(OBJ)/TigrFile.o \
		$(OBJ)/TigrConfigFile.o \
		$(OBJ)/TigrStacktrace.o \
		$(OBJ)/HMMbuilder.o \
		$(OBJ)/HMMreader.o \
		$(OBJ)/Symbol.o \
		$(LIBS)
#--------------------------------------------------------

#--------------------------------------------------------
$(OBJ)/viterbi-trainer.o:\
		viterbi-trainer.C
	$(CC) $(CFLAGS) -o $(OBJ)/viterbi-trainer.o -c \
		viterbi-trainer.C
#---------------------------------------------------------
viterbi-trainer: \
		$(OBJ)/viterbi-trainer.o \
		$(OBJ)/TigrGffReader.o \
		$(OBJ)/HigherOrderAlphabet.o \
		$(OBJ)/NthOrderStringIterator.o \
		$(OBJ)/FastForward.o \
		$(OBJ)/TigrGffExon.o \
		$(OBJ)/TigrGffFeature.o \
		$(OBJ)/TigrGffTranscript.o \
		$(OBJ)/TigrExceptions.o \
		$(OBJ)/TigrBitSet.o \
		$(OBJ)/TigrCommandLine.o \
		$(OBJ)/TigrString.o \
		$(OBJ)/TigrFastaReader.o \
		$(OBJ)/DnaAlphabet.o \
		$(OBJ)/TigrRegex.o \
		$(OBJ)/Alphabet.o \
		$(OBJ)/HiddenMarkovModel.o \
		$(OBJ)/HMMGraph.o \
		$(OBJ)/Sequence.o \
		$(OBJ)/TigrProgress.o \
		$(OBJ)/FastViterbi.o \
		$(OBJ)/TigrRandom.o \
		$(OBJ)/TigrStrTokenizer.o \
		$(OBJ)/TigrFile.o \
		$(OBJ)/TigrConfigFile.o \
		$(OBJ)/TigrStacktrace.o \
		$(OBJ)/HMMbuilder.o \
		$(OBJ)/HMMreader.o \
		$(OBJ)/Symbol.o
	$(CC) $(LDFLAGS) -o viterbi-trainer \
		$(OBJ)/viterbi-trainer.o \
		$(OBJ)/TigrGffReader.o \
		$(OBJ)/HigherOrderAlphabet.o \
		$(OBJ)/NthOrderStringIterator.o \
		$(OBJ)/FastForward.o \
		$(OBJ)/TigrGffExon.o \
		$(OBJ)/TigrGffFeature.o \
		$(OBJ)/TigrGffTranscript.o \
		$(OBJ)/TigrExceptions.o \
		$(OBJ)/TigrBitSet.o \
		$(OBJ)/TigrCommandLine.o \
		$(OBJ)/TigrString.o \
		$(OBJ)/TigrFastaReader.o \
		$(OBJ)/DnaAlphabet.o \
		$(OBJ)/TigrRegex.o \
		$(OBJ)/Alphabet.o \
		$(OBJ)/HiddenMarkovModel.o \
		$(OBJ)/HMMGraph.o \
		$(OBJ)/Sequence.o \
		$(OBJ)/TigrProgress.o \
		$(OBJ)/FastViterbi.o \
		$(OBJ)/TigrRandom.o \
		$(OBJ)/TigrStrTokenizer.o \
		$(OBJ)/TigrFile.o \
		$(OBJ)/TigrConfigFile.o \
		$(OBJ)/TigrStacktrace.o \
		$(OBJ)/HMMbuilder.o \
		$(OBJ)/HMMreader.o \
		$(OBJ)/Symbol.o \
		$(LIBS)

#---------------------------------------------
