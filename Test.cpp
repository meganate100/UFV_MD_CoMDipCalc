#include <iostream>
#include <string>
#include <vector>
#include <fstream>
#include <sstream>
//TestBestTheRest
int main() {
struct Particle {
    float xCoM, yCoM, zCoM, xDip, yDip, zDip; //This will be were we store the final data for 1 particle. Right now I'm just using it as a general vector.
};

std::vector<Particle> Part;
std::ifstream file("grotestcopy.txt");
std::string line;
std::getline(file, line);//We will actually need to consider these lines eventually, but thats another problem for later.
std::getline(file, line);

    std::string DontCare;
    std::string KindaCare;
    float xCoM;
    float yCoM;
    float zCoM;
    float xDip;
    float yDip;
    float zDip;
    float IDK;

for(int i = 0; i < 4; ++i){
std::getline(file, line);
std::istringstream iss(line);
std::cout << std::endl << line;
    iss >> DontCare >> KindaCare >> xCoM >> yCoM >> zCoM >> xDip >> yDip >> zDip >> IDK;//Right now I'm just demoing that this is how we could parse the lines, not actually assigning them to the proper variables
    Part.push_back({xCoM, yCoM, zCoM, xDip, yDip, zDip});
}
std::cout << std::endl << Part.at(0).yCoM;
    return 0;

}

