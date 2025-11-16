#include <iostream>
#include <string>
#include <vector>
#include <fstream>
#include <sstream>
#include <cmath>
//TestBestTheRest
float DotProduct(float x1, float y1, float z1, float x2, float y2, float z2){//Function that will take a string which describes the particle type (eg: "O" for oxygen, "H" for Hydrogen) and translate it to the appropriate mass
float result = ((x1*x2) + (y1*y2) + (z1*z2));
return result;
};

int main() {

struct Particle {
    float timestamp, xCoM, yCoM, zCoM, CoMDist, xDip, yDip, zDip; //This will be were we store the final data for 1 particle. Right now I'm just using it as a general vector.
};

    std::vector<Particle> Part;
    std::ifstream file("droplet_clean.gro");
    std::ofstream TestOutput("TestOutput.txt");
    std::string line;
    std::string DontCare;
    std::string PartType;
    std::istringstream streamMcgee(line);
    int PartNum;
    const float Charge = 0.6732053036584803;
    float xCoM, yCoM, zCoM;
    float CoMDist;
    float H1xCoM, H1yCoM, H1zCoM;
    float H2xCoM, H2yCoM, H2zCoM;
    float VxCoM, VyCoM, VzCoM;
    float xDip, yDip, zDip;
    float timestamp;
    int j = -1;
    int k = 0;

while(std::getline(file, line)){
    streamMcgee.seekg(0, std::stringstream::beg);
    streamMcgee.str(line);
    streamMcgee >> DontCare;

    if(DontCare == "Generated"){//This checks to see if the next line is starting a new timestamp.
        j++;
        streamMcgee.seekg(31, std::stringstream::beg);//Skips to where the timestamp will be in the line
        streamMcgee >> timestamp;
        TestOutput << std::endl << timestamp << std::endl;     
        std::getline(file, line);
        streamMcgee.str(line);
        streamMcgee.seekg(0, std::stringstream::beg);
        streamMcgee >> PartNum;
        TestOutput << PartNum << std::endl;

        for (int i = 1; i < (PartNum/4) +1; i++) {
            TestOutput << std::endl << i;
            std::getline(file, line);//Get Info on the Oxygen
            streamMcgee.str(line);
            TestOutput << line << std::endl;
            streamMcgee.seekg(20, std::stringstream::beg);
            streamMcgee >> xCoM >> yCoM >> zCoM;

            CoMDist = sqrt(DotProduct(xCoM, yCoM, zCoM, xCoM, yCoM, zCoM));//Find radial distance from the CoM of the system.

            std::getline(file, line);//Get Info on the 1st Hydrogen
            streamMcgee.str(line);
            TestOutput << line << std::endl;
            streamMcgee.seekg(20, std::stringstream::beg);
            streamMcgee >> H1xCoM >> H1yCoM >> H1zCoM;

            std::getline(file, line);//Get Info on the 2nd Hydrogen
            streamMcgee.str(line);
            TestOutput << line << std::endl;
            streamMcgee.seekg(20, std::stringstream::beg);
            streamMcgee >> H2xCoM >> H2yCoM >> H2zCoM;

            std::getline(file, line);//Get Info on the virtual site
            streamMcgee.str(line);
            TestOutput << line << std::endl;
            streamMcgee.seekg(20, std::stringstream::beg);
            streamMcgee >> VxCoM >> VyCoM >> VzCoM;

            xDip = (Charge)*(H1xCoM + H2xCoM- 2*VxCoM);//Find Net Dipole moment of Molecule
            yDip = (Charge)*(H1yCoM + H2yCoM- 2*VyCoM);
            zDip = (Charge)*(H1zCoM + H2zCoM- 2*VzCoM);

            Part.push_back({timestamp, xCoM, yCoM, zCoM, CoMDist, xDip, yDip, zDip});//Input all relevant info into a particle structure
            TestOutput << "Timestamp: " << Part.at(((j)*(PartNum/4))+i-1).timestamp << "  xCoM: " << Part.at(((j)*(PartNum/4) )+i-1).xCoM 
            << "  yCoM: " << Part.at(((j)*(PartNum/4) )+i-1).yCoM << "  zCoM: " << Part.at(((j)*(PartNum/4) )+i-1).zCoM 
            << "  Radial Distance: " << Part.at(((j)*(PartNum/4) )+i-1).CoMDist << "  Net xDipole: " 
            << Part.at(((j)*(PartNum/4) )+i-1).xDip << "  Net yDipole: " << Part.at(((j)*(PartNum/4) )+i-1).yDip 
            << "  Net zDipole: " << Part.at(((j)*(PartNum/4) )+i-1).zDip << std::endl;
        }
        
    } else{
        k++;
        std::cout << "The else has been activated " << k <<std::endl;
    }
}
TestOutput << "The Program might have worked" << std::endl;
return 0;
}
