LIBDIR=lib
INCDIR=src
MODELLIB=model_sm
CXXFLAGS= -std=c++1z -O -I$(INCDIR) -I.
LIBFLAGS= -L$(LIBDIR) -l$(MODELLIB)

all: generate_vegas generate_flat convert_to_ps compute_metrics_from_weights compute_metrics_from_likelihoods

generate_vegas: generate_vegas.cpp CPPProcess.o common_elements.hh $(LIBDIR)/lib$(MODELLIB).a
	$(CXX) $(CXXFLAGS) -o $@ generate_vegas.cpp CPPProcess.o $(LIBFLAGS)

generate_flat: generate_flat.cpp CPPProcess.o common_elements.hh $(LIBDIR)/lib$(MODELLIB).a
	$(CXX) $(CXXFLAGS) -o $@ generate_flat.cpp CPPProcess.o $(LIBFLAGS)

convert_to_ps: convert_to_ps.cpp CPPProcess.o common_elements.hh $(LIBDIR)/lib$(MODELLIB).a
	$(CXX) $(CXXFLAGS) -o $@ convert_to_ps.cpp CPPProcess.o $(LIBFLAGS)

compute_metrics_from_weights: compute_metrics_from_weights.cpp CPPProcess.o common_elements.hh $(LIBDIR)/lib$(MODELLIB).a
	$(CXX) $(CXXFLAGS) -o $@ compute_metrics_from_weights.cpp CPPProcess.o $(LIBFLAGS)

compute_metrics_from_likelihoods: compute_metrics_from_likelihoods.cpp CPPProcess.o common_elements.hh $(LIBDIR)/lib$(MODELLIB).a
	$(CXX) $(CXXFLAGS) -o $@ compute_metrics_from_likelihoods.cpp CPPProcess.o $(LIBFLAGS)

$(LIBDIR)/lib$(MODELLIB).a:
	rm -rf lib
	mkdir lib
	cd $(INCDIR) && make

.PHONY: clean

clean:
	rm -f generate_vegas generate_flat convert_to_ps compute_metrics_from_weights compute_metrics_from_likelihoods
	rm -f CPPProcess.o
	rm -rf events
	rm -rf lib
	cd $(INCDIR) && make clean