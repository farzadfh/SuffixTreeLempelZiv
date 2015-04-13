// Gives the suffix tree and Lempel-Ziv factorization of the input string.
// Farzad Farnoud, September 2014
// Uses code from Mark Nelson

#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <string.h>
#include <cassert>
#include <string>
#include <fstream>
#include <algorithm>
#include <vector>

using namespace std;

class Node;

class Edge {
public:
	int first_char_index;
	int last_char_index;
	Node *start_node;
	Node *end_node;
	int id;
	Edge(){
		start_node = NULL;
		end_node = NULL;
		first_char_index = -1;
		last_char_index = -1;
		id = -1;
	}
};

class Node {
public :
	Edge *parent_edge;
	int n_childern; // number of children
	vector<Edge*> child_edges;
	Node *suffix_node;
	int char_depth;
	int id;
	Node(int alphabet_size){
		parent_edge = NULL;
		n_childern = 0;
		child_edges.resize(alphabet_size);
		suffix_node = NULL;
		char_depth = 0;
		id = -1;
	}
};

class SuffixTree {
public:
	string &T; // the sequence for which the suffix tree is constructed
	int N; // length of T
	int M; // size of the alphabet of T
	vector<Edge> Edges;
	Edge * new_edge;
	vector<Node> Nodes;
	Node * new_node;
	vector <int> LSP;
	vector <int> LZF;
	SuffixTree(string &t):T(t){
		N = T.length();
		set_M__unique_last_char(); // ensures uniqueness of last character and sets M
		Edges.resize(2*N);
		new_edge = &Edges[0];
		Nodes.resize(2*N, Node(M));
		new_node = &Nodes[0];
		LSP.resize(N);
	}
	void print_edges(ostream &);
	void build_stree_lzf();
private:
	Edge * find_matching_edge(Node *, char);
	void new_leaf(Node *, int);
	Node * split_edges(Edge *, int);
	void update_suffixless_node(Node*& sfxless, Node* candidate_sfx, Node* new_sfxless);
	void set_M__unique_last_char();
};



// count the number of unique characters in the input string to determine the maximum degree of each node and ensure uniqueness of last character
void SuffixTree::set_M__unique_last_char(){
	// checking last character to be unique
	cout<<"\nChecking uniqueness of last character..."<<flush;
	vector <int> has_char (UCHAR_MAX+1,0);
	bool unique_last_symbol=false;
	int seen_before;
	M = 0;
	for (int i = 0; i<N; ++i) {
		seen_before = has_char[T[i]];
		has_char[T[i]] = 1;
		M += 1-seen_before;
	}
	unique_last_symbol = !seen_before;
	// if the last symbol is not unique try to add $, if it is not used
	if (unique_last_symbol){
		cout<<"Last character unique!\n\n"<<flush;
	} else {
		if (!has_char['$']){
			T.push_back('$');
			++N;
			++M;
			cout<<"Last character not unique. Adding '$'.\n\n"<<flush;
		} else {
			throw std::invalid_argument("Last character not unique. Cannot add '$'.\nEnd string with unique symbol or ensure it does not contain '$'.\n");

		}
	}
}

Edge * SuffixTree::find_matching_edge(Node * ndp, char c){
	for (int i = 0; i < ndp->n_childern; ++i)
		if (T[ndp->child_edges[i]->first_char_index + 1] == c){
			return ndp->child_edges[i];
			break;
		}
		return NULL;
}



// show the grahp
void SuffixTree::print_edges(ostream& of){
	of<< " Start  End  Suf  First Last  String\n";
	for (int i = 0; ; ++i){
		Edge *e = &Edges[0] + i;
		if (! e->start_node )
			break;
		of << setw( 5 ) << e->start_node->id << " "
			<< setw( 5 ) << e->end_node->id << " "
			<< setw( 3 ) << (e->end_node->suffix_node ? e->end_node->suffix_node->id : -1) << " "
			<< setw( 5 ) << e->first_char_index << " "
			<< setw( 6 ) << e->last_char_index << "  ";
		for ( int l = e->first_char_index+1 ;
			l <= e->last_char_index;
			l++ )
			of << T[ l ];
		of << "\n";
		of.flush();
	}
}

Node * SuffixTree::split_edges(Edge * edp, int split_point){ // split_point is the number characters that are given to the new edge

	// set up new edge
	new_edge->first_char_index = edp->first_char_index;
	new_edge->start_node = edp->start_node;
	new_edge->end_node = new_node;
	new_edge->last_char_index = edp->first_char_index + split_point;
	new_edge->id = new_edge - &Edges[0];

	// update old parent
	for (int l = 0; l < edp->start_node->n_childern; ++l)
		if (edp->start_node->child_edges[l] == edp)
			edp->start_node->child_edges[l] = new_edge;

	// set up new node
	new_node->char_depth = new_edge->start_node->char_depth + split_point;
	new_node->child_edges[0] = edp;
	new_node->n_childern = 1;
	new_node->parent_edge = new_edge;
	new_node->id = new_node - &Nodes[0];

	// update old edge
	edp->first_char_index = edp->first_char_index + split_point;
	edp->start_node = new_node;

	new_edge++;
	return new_node++;
}

void SuffixTree::new_leaf(Node * ndp, int k){
	new_edge->start_node = ndp;
	new_edge->end_node = new_node;
	new_edge->first_char_index = k;
	new_edge->last_char_index = N-1;
	new_edge->id = new_edge - &Edges[0];
	new_node->char_depth = ndp->char_depth + N-1 - k;
	new_node->n_childern = 0;
	new_node->parent_edge = new_edge;
	new_node->id = new_node - &Nodes[0];

	ndp->child_edges[ndp->n_childern++] = new_edge;
	new_node++;
	new_edge++;
}


void SuffixTree::build_stree_lzf()
{
	cout<<"Initializing suffix tree..."<<flush;
	Node * current_node;
	int current_offset;
	Edge * current_edge;

	Node * root = new_node; // set up root
	root->suffix_node = root; // suffix node of root is itself
	root->id = 0;
	new_node++;
	current_node = root;
	current_offset = 0;
	int i = 0;
	int j = -1;
	int k = -1;
	Edge * edge_ix = NULL;
	Node * suffixless_node = NULL;
	int edge_length;
	cout<< "done!\n\n" << flush;
	cout<< "Constructing suffix tree of sequence of length "<< N << "..." << flush;
	while (i<N){
		if (current_offset && current_edge->first_char_index + current_offset == current_edge->last_char_index) {
			// at node
			current_node = current_edge->end_node;
			current_offset = 0;
		}
		if (current_offset == 0){
			// explicit node
			update_suffixless_node(suffixless_node, current_node, NULL);
			edge_ix = find_matching_edge(current_node,T[k+1]);
			if ( edge_ix == NULL) {
				// no match found, create new leaf
				new_leaf(current_node,k);
				LSP[j+1] = k-j;
				k = ++j; i += (i==j);
				if (current_node != root){
					if (current_node->suffix_node)
						current_node = current_node->suffix_node;
					else
						current_node = current_node->parent_edge->start_node->suffix_node;
					k = j + current_node->char_depth;
				}
			} else {
				// match found, count and skip
				edge_length = edge_ix->last_char_index - edge_ix->first_char_index;
				if (k == i-1){
					current_edge = edge_ix;
					current_offset = 1;
					k = i++;
				} else if(k + edge_length < i) {
					// skip whole edge
					k = k + edge_length;
					current_node = edge_ix->end_node;
				} else {
					// skip part of edge
					current_edge = edge_ix;
					current_offset = i - 1 - k;
					k = i - 1;
				}
			}
		} else {
			// implicit node
			if (T[k+1] == T[current_edge->first_char_index + current_offset+1]){
				// next character matches
				k = i++;
				current_offset++;
			} else {
				// no match, create new node and split_edges edge
				current_node = split_edges(current_edge,current_offset);
				current_offset = 0;
				update_suffixless_node(suffixless_node, current_node, current_node);
			}
		}
	}
	cout<<"done!\n"<<flush;
	cout<<"Suffix tree has " << new_node - &Nodes[0] << " nodes.\n\n" << flush;
	cout<<"Computing LZF..."<<flush;
	for (unsigned long j=0;j<N;){
		LZF.push_back((LSP[j]? LSP[j]:1));
		j += (LSP[j]? LSP[j]:1);
	}
	cout<<"done!\n\n"<<flush;
}

void SuffixTree::update_suffixless_node(Node*& sfxless, Node* candidate_sfx, Node* new_sfxless){
	if (new_sfxless) {
		if (sfxless){
			assert(sfxless->char_depth-1 == candidate_sfx->char_depth);
			sfxless->suffix_node = candidate_sfx;
		}
		sfxless = new_sfxless;
	} else {
		if (sfxless && sfxless->char_depth-1 == candidate_sfx->char_depth){
			sfxless->suffix_node = candidate_sfx;
			sfxless = NULL;
		}
	}
}



void get_sequence(ifstream & InFile, int start_line, int end_line, string &t){
	string str;
	t.assign("");
	int line_num=1;
	int lines_read = 0;
	while (line_num<start_line && getline(InFile, str)) ++line_num; // ignoring lines 1..start_line
	if (end_line){ // if an end line number is given (i.e., it isn't zero), read only until (inclusive) that line
		cout<<"\nReading input file lines " << start_line << " through " << end_line << "..." <<flush;	
		while ( line_num<=end_line && getline(InFile, str)){
			transform(str.begin(),str.end(),str.begin(),::tolower); // convert to lower case
			t += str;
			++line_num;
			++lines_read;
		}
	} else {
		cout<<"\nReading input file lines " << start_line << " through end of file" << "..." <<flush;	
		while (getline(InFile, str)){
			transform(str.begin(),str.end(),str.begin(),::tolower); // convert to lower case
			t += str;
			++line_num;
			++lines_read;
		}
	}
	cout<<lines_read<<" lines read!\n"<<flush;
	InFile.close();
	int n = t.length();
	if (n==0) {
		throw std::invalid_argument("Sequence Empty.\n");
	}
}



int main(int argc, char ** argv) {
	ifstream InFile;
	ofstream SuffixTreeFile, LZFFile;
	int code=0;


	//set default values
	int start_line = 1;
	int end_line = 0;
	string input_sequence_file;
	string stree_file;
	string LZF_file;

	string t;

	int input_mode = 0; // 0: no input given. demo with default. 1: input filename given. 2: input sequence given
	bool isInputTyped = false;
	bool isInputFile = false;
	bool isLSPRequested = false;
	bool isStreeRequested = false;
	bool isLZFRequested = false;
	

	try {
	for (int i = 1; i < argc; ++i) {
		std::string arg = argv[i];
		if ((arg == "-h") || (arg == "--help")){
			cout<<"\nNAME\n";
			cout<<"\tSTree_LZF\n";
			cout<<"\tGives the suffix tree and Lempel-Ziv factorization of the input string\n";
			cout<<"USAGE 1: \n";
			cout<<"\tSTree_LZF -s input-sequence\n" << flush;
			cout<<"USAGE 2: \n";
			cout<<"\tSTree_LZF -f input-file [-b start-line] [-e end-line] [-t suffix-tree-output-file]\n\t[-l LZ-factorization-output-file [-p]] \n" << flush;
			cout<<"NOTES:\n";
			cout<<"\t-b and -e indicate the line numbers within which the input file is processed.\n";
			cout<<"\tIf -p is used, a list is given where the i'th number is the length  of longest \n\t\tprefix starting at i that also starts at j<i.\n";
			cout<<"\tIf -s is used, all outputs are printed on screen.\n"; 
			return 0;
		} else if ((arg == "-f") || (arg == "--input_file")) {
			if ( i + 1 < argc) {
				isInputFile = true;
				input_sequence_file.assign(argv[++i]);
			} else {
				throw std::invalid_argument("No argument given for the option --input_file.\n");
			}
		} else if ((arg == "-s") || (arg == "--input_seq")) { // For small inputs. All outputs are printed on screen.
			if ( i + 1 < argc) {
				isInputTyped = true;
				t.assign(argv[++i]);
			} else {
				throw std::invalid_argument("No argument given for the option --input_seq.\n");
			} 
		} else if ((arg == "-t") || (arg == "--suffix_tree_file")) {
			if ( i + 1 < argc) {
				isStreeRequested = true;
				stree_file.assign(argv[++i]);
			} else {
				throw std::invalid_argument("No argument given for the option --suffix_tree_file.\n");
			} 
		} else if ((arg == "-l") || (arg == "--LZF_file")) {
			if ( i + 1 < argc) {
				isLZFRequested = true;
				LZF_file.assign(argv[++i]);
			} else {
				throw std::invalid_argument("No argument given for the option --LZF_file.\n");
			} 
		} else if ((arg == "-b") || (arg == "--start_line")) {
			if ( i + 1 < argc) {
				start_line = atoi(argv[++i]);
			} else {
				throw std::invalid_argument("No argument given for the option --start_line.\n");
			} 
		} else if ((arg == "-e") || (arg == "--end_line")) {
			if ( i + 1 < argc) {
				end_line = atoi(argv[++i]);
			} else {
				throw std::invalid_argument("No argument given for the option --end_line.\n");
			} 
		} else if ((arg == "-p") || (arg == "--output_LSP")) { // Not described in USAGE
			isLSPRequested = true;
		}
	}


	if (!isInputFile && !isInputTyped){
		throw std::invalid_argument("No input given.\n");
	}

	if (!isLZFRequested && !isStreeRequested && !isInputTyped){
		throw std::invalid_argument("No output requested.\n");
	}

	if(isInputFile){
		InFile.open(input_sequence_file);
		get_sequence(InFile, start_line, end_line, t);
	}


#ifdef NDEBUG
	if (isInputTyped && (isLZFRequested || isStreeRequested)){ // If input is manual, output is only on screen
		cout << "Options -t and -l are not compatible with -s. Ignoring -t and/or -l.\n";
		isLZFRequested = isStreeRequested = false;
	}
#endif
		SuffixTree st(t);
		st.build_stree_lzf();

		// OUTPUT
		if (isLZFRequested){
			LZFFile.open(LZF_file);
			cout<<"Writing LZF to "<<LZF_file<<" ..."<<flush;
			for (vector<int>::const_iterator i = st.LZF.begin(); i != st.LZF.end(); ++i)
				LZFFile << *i << ' ';

			if(isLSPRequested){
				LZFFile<<"\n";
				for (int j=0;j<st.N;++j)
					LZFFile << st.LSP[j]<<" ";
			}
			LZFFile.close();
			cout<<"done!\n\n"<<flush;
		}


		if (isStreeRequested){
			SuffixTreeFile.open(stree_file);
			cout<<"Writing suffix tree to "<<stree_file<<" ..."<<flush;
			st.print_edges(SuffixTreeFile);
			SuffixTreeFile.close();
			cout<<"done!\n\n"<<flush;
		} 

		if (isInputTyped){ // input given as a sequence not file. print output to screen
			cout<<"Suffix tree:\n"<<flush;
			st.print_edges(cout);
			cout<<"done!\n\n"<<flush;
			cout<<"Block lengths in Lempel-Ziv factorization:\n"<<flush;
			for (vector<int>::const_iterator i = st.LZF.begin(); i != st.LZF.end(); ++i)
				cout << *i << ' ';

			cout<<"\n\nThe i'th number is the length of longest prefix starting at i that also starts at j<i:\n";
			for (int j=0;j<st.N;++j)
				cout << st.LSP[j]<<" ";
			cout<<"\ndone!\n\n"<<flush;
		}

	}
	catch( const std::invalid_argument& e ) {
		cout<<e.what();
		code = 1;
	}
	cout<<"Press Enter to exit..."<<flush;
	cin.get();
	cout<<"\n";
	return(code);
}

