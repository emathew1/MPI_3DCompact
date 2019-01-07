#include <iostream>
#include <fstream>
#include <string>
#include <boost/program_options.hpp>

using namespace std;
namespace po = boost::program_options;

int main(int argc, char *argv[]){

    ifstream config_file("test.in");

    po::options_description input("Input Parameters");

    int test1a = -1, test1b = -1;
    int test2 =-1;

    input.add_options()
	("test1.test1a", po::value<int>(&test1a), "test_item_1a")
	("test1.test1b", po::value<int>(&test1b), "test_item_1b")
	("test2.test2",  po::value<int>(&test2), "test_item_2");

    po::variables_map vm;
    po::store(po::parse_config_file(config_file, input), vm);
    po::notify(vm);
	

    if(vm.count("test1.test1a")){
	cout << " test1.test1a = " << test1a <<  endl;
    }else{
	cout << " Didn't find test1.test1a value " << endl;
    }

    if(vm.count("test1.test1b")){
	cout << " test1.test1b = " << test1b <<  endl;
    }else{
	cout << " Didn't find test1.test1b value " << endl;
    }



    if(vm.count("test2.test2")){
	cout << " test2 = " << test2 << endl;
    }else{
	cout << " Didn't find test2 value " << endl;
    }



};
