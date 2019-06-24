#ifndef JuniprJetLoader_H
#define JuniprJetLoader_H
#include <iostream>

#include "JuniprUtils.hh"

using namespace fastjet;
using namespace std;

class JuniprJet
{
public:
  JuniprJet(vector< vector<double> > CSJets_arg,
	    vector< vector<int> >    intermediate_states_arg, 
	    vector<int>              mothers_arg, 
	    vector< vector<int> >    daughters_arg,
	    vector< vector<double> > branchings_arg,
	    int label_arg){
    label               = label_arg;
    CSJets              = CSJets_arg;
    intermediate_states = intermediate_states_arg;
    mothers             = mothers_arg;
    daughters           = daughters_arg;
    branchings          = branchings_arg;
    
    set_multiplicity();
    set_n_branchings();
    set_seed_momentum();
    set_mother_momenta();
    set_daughter_momenta();
    set_mother_id_energy_order();
  }
  ~JuniprJet() {}
        
  // Variables set by constructor:
  
  // CSJets contains the ClusteringSequenceJets. Each element is a 4-momentum (E, theta, phi, mass) 
  vector< vector<double> > CSJets;
  
  //// intermediate_states, mother and daugthers contain indices to the entries in CSJets
  // intermediate_states[i] is a list of the intermediate states at step i in energy order (high to low)
  vector< vector<int> >    intermediate_states;
  // mothers[i] is the index to the decaying mother at time step i
  vector<int>              mothers;
  // daughters[i] is the index to the two daughters at time step i. Most energetic daughter first. 
  vector< vector<int> >    daughters;
  // branchings[i] is the branchings (z, theta, phi, delta) at time step i. 
  vector< vector<double> > branchings;
  
  // mother_id_energy_order[i] is the id of the mother that branches in the list of energy ordered intermediate states at time step i
  vector<int> mother_id_energy_order;
  
  int label; // label for used for classification target value
  
  int multiplicity; // number of final state particles
  int n_branchings; // number of branchings/
  
  // Seed_momentum = sum of all final state momenta
  vector<double>                     seed_momentum;
  // daughter_momenta[i] contains the two 4-momenta (E, theta, phi, mass) of the daughters produced at time step i.
  vector< vector< vector<double> > > daughter_momenta;
  // mother_momenta[i] contains the 4-momenta (E, theta, phi, mass) of the mother that branches at time step i.
  vector< vector<double> >           mother_momenta;
  
  void write_to_json(ostream& outfile, vector < double > vars){
    outfile << "{\n\t" ;
    json_label(outfile);
    outfile << ",\n\t";
    json_vars(outfile,vars);
    outfile << ",\n\t";
    json_multiplicity(outfile);
    outfile << ",\n\t";
    json_n_branchings(outfile);
    outfile << ",\n\t";
    json_seed_momentum(outfile);
    outfile << ",\n\t";
    json_CSJets(outfile);
    outfile << ",\n\t";
    json_intermediate_states(outfile);
    outfile << ",\n\t";
    json_mother_id_energy_order(outfile);
    outfile << ",\n\t";
    json_mothers(outfile);
    outfile << ",\n\t";
    json_daughters(outfile);
    outfile << ",\n\t";
    json_branchings(outfile);
    outfile << ",\n\t";
    json_mother_momenta(outfile);
    outfile << ",\n\t";
    json_daughter_momenta(outfile);
    outfile << "\n}" ;
  }  
  
private:
  void set_multiplicity(){ multiplicity = (CSJets.size()+1)/2; }
  void set_n_branchings(){ n_branchings = (CSJets.size()-1)/2;}
  void set_seed_momentum(){seed_momentum = CSJets[mothers[0]];}
  void set_mother_momenta(){
    for(int i = 0; i< int(mothers.size()); i++){
      mother_momenta.push_back(CSJets[mothers[i]]);
    }
  }
  void set_daughter_momenta(){
    for(int i = 0; i< int(daughters.size()); i++){
      vector<double> daughter1 = CSJets[daughters[i][0]];
      vector<double> daughter2 = CSJets[daughters[i][1]];
      
      vector< vector<double> > two_daughters;
      two_daughters.push_back(daughter1);
      two_daughters.push_back(daughter2);
      daughter_momenta.push_back(two_daughters);
    }
  }
  void set_mother_id_energy_order(){
    for (int i = 0; i< int(mothers.size()); i++){
      for (int j = 0; j< int(intermediate_states[i].size()); j++){
	if(mothers[i] == intermediate_states[i][j]){
	  mother_id_energy_order.push_back(j);
	  break;
	}
      }
    }
    if(mothers.size() !=mother_id_energy_order.size()){
      cout << "Error in set_mother_id_energy_order(): one or more mother_ids not found in intermediate_states" << endl;
    }
  }
    
  // Printing functions in json format
  void json_momentum(ostream& outfile, vector<double> momentum){
    outfile << "[ " 
            << momentum[0] << " , "
            << momentum[1] << " , "
            << momentum[2] << " , "
            << momentum[3] << " ]";
  }
  void json_vector_int(ostream& outfile, vector<int> vector_int){
    outfile << "[ ";
    
    for(int i = 0; i<int(vector_int.size()); i++){
      outfile << vector_int[i];
      if(i!=int(vector_int.size()-1)){
	outfile << " , ";
      }      
    }
    outfile << " ]";
  }
  void json_vector_double(ostream& outfile, vector<double> vector_double){
    outfile << "[ ";
    
    for(int i = 0; i<int(vector_double.size()); i++){
      outfile << vector_double[i];
      if(i!=int(vector_double.size()-1)){
	outfile << " , ";
      }      
    }
    outfile << " ]";
  }
  void json_label(ostream& outfile){ outfile << "\"label\" : " << label ;}
  void json_vars(ostream& outfile, vector<double> vars){ 
    outfile << "\"vars\" : ";
    json_vector_double(outfile, vars);
  }
  void json_multiplicity(ostream& outfile){ outfile << "\"multiplicity\" : " << multiplicity ;}
  void json_n_branchings(ostream& outfile){  outfile << "\"n_branchings\" : " << n_branchings ; }
  void json_seed_momentum(ostream& outfile){ 
    outfile << "\"seed_momentum\" :" ;
    json_momentum(outfile, seed_momentum);
  }
  void json_CSJets(ostream& outfile){
    outfile << "\"CSJets\" : [\n\t\t";
    for (int i = 0; i<int(CSJets.size()); i++){
      json_momentum(outfile, CSJets[i]);
      if(i!=int(CSJets.size()-1)){
	outfile << ",\n\t\t";
      }
    }
    outfile << "]";
  }
  void json_intermediate_states(ostream& outfile){
    outfile << "\"CS_ID_intermediate_states\": [\n\t\t";
    for (int i = 0; i< int(intermediate_states.size()); i++){
      json_vector_int(outfile, intermediate_states[i]);
      if(i!=int(intermediate_states.size()-1)){
	outfile << ",\n\t\t";
      }
    }
    outfile << "]";
  }
  void json_mothers(ostream& outfile){
    outfile << "\"CS_ID_mothers\": ";
    json_vector_int(outfile, mothers);
  }
  void json_mother_id_energy_order(ostream& outfile){
    outfile << "\"mothers_id_energy_order\": ";
    json_vector_int(outfile, mother_id_energy_order);
  }
  void json_daughters(ostream& outfile){
    outfile << "\"CS_ID_daughters\": [\n\t\t";
    for (int i = 0; i< int(daughters.size()); i++){
      json_vector_int(outfile, daughters[i]);
      if(i!=int(daughters.size()-1)){
	outfile << ",\n\t\t";
      }
    }
    outfile << "]";
  }
  void json_branchings(ostream& outfile){
    outfile << "\"branchings\" : [\n\t\t";
    for (int i = 0; i<int(branchings.size()); i++){
      json_vector_double(outfile, branchings[i]);
      if(i!=int(branchings.size()-1)){
	outfile << ",\n\t\t";
      }
    }
    outfile << "]";
  }
  void json_mother_momenta(ostream& outfile){
    outfile << "\"mother_momenta\" : [\n\t\t";
    for (int i = 0; i<int(mother_momenta.size()); i++){
      json_momentum(outfile, mother_momenta[i]);
      if(i!=int(mother_momenta.size()-1)){
	outfile << ",\n\t\t";
      }
    }
    outfile << "]";
  }
  void json_daughter_momenta(ostream& outfile){
    outfile << "\"daughter_momenta\" : [\n\t\t";
    for (int i = 0; i<int(daughter_momenta.size()); i++){
      outfile << "[";
      json_momentum(outfile, daughter_momenta[i][0]);
      outfile << " , ";
      json_momentum(outfile, daughter_momenta[i][1]);
      outfile << " ]";

      if(i!=int(daughter_momenta.size()-1)){
	outfile << ",\n\t\t";
      }
    }
    outfile << "]";
  }


};

class JuniprJetLoader {
public:
  JuniprJetLoader(ClusterSequence cs, int label);
  ~JuniprJetLoader();
  vector<double> get_branching(vector<PseudoJet> const &PJ_CSJets, int hist_id, vector<int> current_daughters);
  void remove_hist_id(vector<int> &current_intermediate_state, int hist_id);
  void add_hist_id(vector<int> &current_intermediate_state, int hist_id, vector< vector<double> > const &CSJets);
  void cluster_sequence_to_JuniprJet(ClusterSequence cs, int label);
  void write_juniprjet(ostream& outfile,vector < double > vars);
  JuniprJet *juniprjet;
};
#endif
