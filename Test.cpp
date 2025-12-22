#include <iostream>
#include <string>
#include <vector>
#include <fstream>
#include <sstream>
#include <cmath>
// Structs
struct Particle {
    float timestamp, x, y, z, CoMDist, xDip, yDip, zDip; //This will be where we store the final data for 1 particle. Right now I'm just using it as a general vector.
};

struct Shell {
    float timestamp, NetXDip, NetYDip, NetZDip; //This will be where we store the timestamp and net dipol moment for a given shell.
};

//Test/Temporary Variables
std::vector<std::vector<Shell>> TwoDArray;
std::vector<Shell> TwoDArrayInput;
float xBoxSize;
float yBoxSize;
float zBoxSize;
bool ShellFlag = true;
bool PartNumFlag = true;


// Variable Declaration
std::vector<Particle> Part;
std::string line;
std::string DontCare;
std::string PartType;
std::istringstream stream(line);
int PartNum;
int MolNum;
const float Charge = 0.6732053036584803;
float radial;
int ShellNum;
float xCoM = 0, yCoM = 0, zCoM = 0;
float x, y, z;
float CoMDist = 0;
float H1xCoM, H1yCoM, H1zCoM;
float H2xCoM, H2yCoM, H2zCoM;
float VxCoM, VyCoM, VzCoM;
float xDip, yDip, zDip;
float timestamp;
int j = -1;
int k = 0;

//Function that will take a string which describes the particle type (eg: "O" for oxygen, "H" for Hydrogen) and translate it to the appropriate mass
float DotProduct(float x1, float y1, float z1, float x2, float y2, float z2)  {   
    float result = ((x1*x2) + (y1*y2) + (z1*z2));
    return result;
}

void CoMAdjustment(int j) {
    std::cout << "CoM function has been entered" << std::endl;
    for (int i = 1; i < MolNum +1; i++) {
        x = Part.at(i-1).x - xCoM;//Get the position vectors of the particles to be relative to the CoM of the system.
        y = Part.at(i-1).y - yCoM;
        z = Part.at(i-1).z - zCoM;

        CoMDist = sqrt(DotProduct(x, y, z, x, y, z));//Find radial distance from the CoM of the system.

        Part.at(i-1).x = x;//Re-enters the position values relative to the CoM of the system for all the particles
        Part.at(i-1).y = y;                
        Part.at(i-1).z = k;
        Part.at(i-1).CoMDist = CoMDist;

        for(int n = 1; n <= ShellNum; n++) {//This creates the net dipole moment of a radial shell of a given timestamp.
            if((CoMDist/radial) <= n){
                //std::cout << "The current Net X Dipole moment of Shell #"<< n <<" is: " << TwoDArrayInput.at(n-1).NetXDip <<
                //" and the next particle's x Dipole moment is: " << Part.at(((j)*MolNum)+i-1).xDip << std::endl;

                TwoDArrayInput.at(n-1).NetXDip = TwoDArrayInput.at(n-1).NetXDip + Part.at(i-1).xDip;
                TwoDArrayInput.at(n-1).NetYDip = TwoDArrayInput.at(n-1).NetYDip + Part.at(i-1).yDip;
                TwoDArrayInput.at(n-1).NetZDip = TwoDArrayInput.at(n-1).NetZDip + Part.at(i-1).zDip;
                break;
            }
        } 
    }

    for (int n = 1; n <= ShellNum; n++) {//Adds the timestamp to the shells
        TwoDArrayInput.at(n-1).timestamp = Part.at(0).timestamp;
    }

    TwoDArray.push_back({TwoDArrayInput});


    for(int n = 1; n <= ShellNum; n++) {
        TwoDArrayInput.at(n-1).timestamp = 0;
        TwoDArrayInput.at(n-1).NetXDip = 0;
        TwoDArrayInput.at(n-1).NetYDip = 0;
        TwoDArrayInput.at(n-1).NetZDip = 0;
    }
    std::cout << "Timestamp completed for CoM adjustment function: " << TwoDArray.at(j).at(0).timestamp << std::endl;
    return;
}

int GetFiles()  {
    return 0;
}

int main(int argc, char* argv[])  {
    if (argc != 4)  {
        std::cout << "Please pass the input and output filename arguments respectively, and then the radius of the shell size for comparison. \nExample: test <input-filename> <out-filename> <shell-size-radius>" << std::endl;
        return 0;      
    }

    radial = atof(argv[3]);
    std::ifstream file(argv[1]);
    std::ofstream TestOutput(argv[2]);  

    while(std::getline(file, line)){
        stream.seekg(0, std::stringstream::beg);
        stream.str(line);
        stream >> DontCare;

        if(DontCare == "Generated"){//This checks to see if the next line is starting a new timestamp.
            j++;
            stream.seekg(31, std::stringstream::beg);//Skips to where the timestamp will be in the line
            stream >> timestamp;
            TestOutput << std::endl << timestamp << std::endl;     
            std::getline(file, line);
            stream.str(line);
            stream.seekg(0, std::stringstream::beg);
            stream >> PartNum;
            TestOutput << PartNum << std::endl;
            MolNum = (PartNum/4);

            if(PartNumFlag) {
                Part = std::vector<Particle>(MolNum, {0,0,0,0,0,0,0,0});
                PartNumFlag = false;
                std::cout << Part.size();;
            }

            for (int i = 1; i < MolNum +1; i++) {
                TestOutput << std::endl << i;
                std::getline(file, line);//Get Info on the Oxygen
                stream.str(line);
                TestOutput << line << std::endl;
                stream.seekg(20, std::stringstream::beg);
                stream >> x >> y>> z;
                xCoM = xCoM + x; //Builds the center of mass vector for the configuration
                yCoM = yCoM + y;
                zCoM = zCoM + z;

                std::getline(file, line);//Get Info on the 1st Hydrogen
                stream.str(line);
                TestOutput << line << std::endl;
                stream.seekg(20, std::stringstream::beg);
                stream >> H1xCoM >> H1yCoM >> H1zCoM;

                std::getline(file, line);//Get Info on the 2nd Hydrogen
                stream.str(line);
                TestOutput << line << std::endl;
                stream.seekg(20, std::stringstream::beg);
                stream >> H2xCoM >> H2yCoM >> H2zCoM;

                std::getline(file, line);//Get Info on the virtual site
                stream.str(line);
                TestOutput << line << std::endl;
                stream.seekg(20, std::stringstream::beg);
                stream >> VxCoM >> VyCoM >> VzCoM;

                xDip = (Charge)*(H1xCoM + H2xCoM- 2*VxCoM);//Find Net Dipole moment of Molecule
                yDip = (Charge)*(H1yCoM + H2yCoM- 2*VyCoM);
                zDip = (Charge)*(H1zCoM + H2zCoM- 2*VzCoM);


                Part.at(i-1) = {timestamp, x, y, z, CoMDist, xDip, yDip, zDip};//Input all relevant info into a particle structure
                TestOutput << "Timestamp: " << Part.at(i-1).timestamp << "  x: " << Part.at(i-1).x 
                << "  y: " << Part.at(i-1).y << "  z: " << Part.at(i-1).z << "  Net xDipole: " 
                << Part.at(i-1).xDip << "  Net yDipole: " << Part.at(i-1).yDip 
                << "  Net zDipole: " << Part.at(i-1).zDip << std::endl;
            }   
        } 
        else    {
            k++;

            if(ShellFlag) {//Gets the box size info for the data, which is used for finding how many shells there should be
                std::cout << line;
                stream.str(line);
                stream >> xBoxSize >> yBoxSize >> zBoxSize;
                std::cout << "The box size is: " << xBoxSize << "  " << yBoxSize << "  " << zBoxSize << std::endl; 
                ShellNum = std::ceil((sqrt(DotProduct(xBoxSize, yBoxSize, zBoxSize, xBoxSize, yBoxSize, zBoxSize)))/radial);
                TestOutput << "The Number of Shells is: " << ShellNum;
                TwoDArrayInput = std::vector<Shell>(ShellNum, {0,0,0,0});
               
                ShellFlag = false;
            }
            std::cout << "First pass for timestamp t=" << timestamp << " has completed" <<std::endl;
            xCoM = xCoM/MolNum;//Finalize the real center of mass vector based on all the particles
            yCoM = yCoM/MolNum;
            zCoM = zCoM/MolNum;
            std::cout << "The center of mass for the configuration at t=" << timestamp << " is xCoM: " << xCoM << " yCoM: " << yCoM << " zCoM: " << zCoM << std::endl;
            CoMAdjustment(j);
            xCoM = 0;//Reset CoM variables for the system
            yCoM = 0;
            zCoM = 0;
        }
    }
TestOutput << "The Program might have worked" << std::endl;
return 0;
}
