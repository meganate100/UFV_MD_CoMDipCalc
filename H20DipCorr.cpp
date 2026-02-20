#include <iostream>
#include <string>
#include <vector>
#include <fstream>
#include <sstream>
#include <cmath>

// Structs
struct Particle {
    float timestamp, x, y, z, distCoM, xDip, yDip, zDip; //Stores the relevant data for an H2O particle
};

struct Shell {
    float timestamp, NetXDip, NetYDip, NetZDip; //Stores net dipole moment for a given shell of particles at a particle timestamp
};

// Variable Declaration

const float CHARGE = 0.6732053036584803; //Charge used for H20 Dipole moment

std::string line;
std::string inputString;
std::istringstream stream(line);
std::vector<Particle> particle; //Array for holding all particles in the system
std::vector<std::vector<Shell>> dataArray; //2D Array for holding all the shells' net dipole moment for a given timestamp.
std::vector<Shell> workingArray; //Array for building up data on all the shells' net dipole moment.

int partNum; //Number of particles (Hydrogen, Oxygen, Virtual Sites)
int molNum; //Number of molecules (H2O)
int numOfShells; //Number of shells in the system
int currentShellNum;

float shellWidth; //Radial width of each shell
float xCoM = 0, yCoM = 0, zCoM = 0; //Center of mass compoents
float x, y, z; //Components of molecule/oxygen position
float distCoM = 0; //Distance of a moleculue from the center of mass
float H1x, H1y, H1z;//Components of 1st hydrogen position
float H2x, H2y, H2z; //Components of 2nd hydrogen position
float Vx, Vy, Vz; //Components of the virtual site position
float xDip, yDip, zDip; //Compoents of the dipole moment of the molecule
float timestamp; 
float dt, t0, tf = 0; //Timestep, Intial timestamp, and Final timestamp

//Function computes the dot product between two vectors using the cartesian components of the vectors as inputs
float DotProduct(float x1, float y1, float z1, float x2, float y2, float z2)  {   
    float result = ((x1*x2) + (y1*y2) + (z1*z2));
    return result;
}

// Get max shell size for the system and adjust array size if new shells are required
void CheckShellSize(float CoM) {
    currentShellNum = std::ceil(CoM/shellWidth);
    if (currentShellNum > numOfShells) {
        numOfShells = currentShellNum;
        for (int j = 0; j < dataArray.size(); j++) {
            for (int n = workingArray.size(); n < numOfShells; n++) {
                dataArray.at(j).push_back({0,0,0,0});
            }   
        }
        for (int n = workingArray.size(); n < numOfShells; n++) {
            workingArray.push_back({0,0,0,0});
        }
    }
}

//Calaculates correlation function for all shells describing the system 
void CorrelationCalulation(std::ofstream& out) {
    t0 = dataArray.at(0).at(0).timestamp; //Find the initial timestamp
    dt = dataArray.at(1).at(0).timestamp - t0; //Finds the difference between timesteps
    tf = dataArray.back().at(0).timestamp; //Find the final timestamp
    int n = dataArray.size();
    int k = 0;
    float cor = 0;
    out << "Timesteps: " << dt << "\nInitial Time: " << t0 << "\nFinal Time: " << tf << "\n";////Outputs general information about the system
    out << "Number of Shells: " << numOfShells << "\nShell Width: " << shellWidth << "\n";
    for(int i = 0; i < dataArray.back().size(); i++) { // Goes through all the shells and calculate the correlation function data for each individual shell
        out << "Shell: " << i << "\n"; //Outputs the shell number
        while (n-k > 0) { //Calculates the correlation function value for all applicable timestamps
            for(int j = 0; j < n-k; j++){
                cor = cor + (DotProduct(dataArray.at(j).at(i).NetXDip, dataArray.at(j).at(i).NetYDip, dataArray.at(j).at(i).NetZDip,
                                        dataArray.at(j+k).at(i).NetXDip, dataArray.at(j+k).at(i).NetYDip, dataArray.at(j+k).at(i).NetZDip))/(n-k);
            }
            out << t0 + k*dt << " " << cor << "\n"; //Outputs the timestamp and the corresponding value of the correlation function
            cor = 0;
            k++;
        }
        k = 0;
    }
}

// Takes the calulated center of mass of the system and readjusts all the particle's positions relative to the center of mass
// Sorts system into their corresponding shells based on their distance from the center of mass
void CoMAdjustment(std::ofstream& test) {
    for (int i = 1; i < molNum +1; i++) {
        x = particle.at(i-1).x - xCoM;//Get the position vectors of the particles to be relative to the CoM of the system
        y = particle.at(i-1).y - yCoM;
        z = particle.at(i-1).z - zCoM;

        distCoM = sqrt(DotProduct(x, y, z, x, y, z));//Find shellWidth distance from the CoM of the system

        CheckShellSize(distCoM);

        particle.at(i-1).x = x;//Re-enters the position values relative to the CoM of the system for all the particles
        particle.at(i-1).y = y;                
        particle.at(i-1).z = z;
        particle.at(i-1).distCoM = distCoM;

        for(int n = 1; n <= numOfShells; n++) {//This creates the net dipole moment of a shellWidth shell of a given timestamp
            if((distCoM/shellWidth) <= n){
                workingArray.at(n-1).NetXDip = workingArray.at(n-1).NetXDip + particle.at(i-1).xDip;
                workingArray.at(n-1).NetYDip = workingArray.at(n-1).NetYDip + particle.at(i-1).yDip;
                workingArray.at(n-1).NetZDip = workingArray.at(n-1).NetZDip + particle.at(i-1).zDip;
                break;
            }
        } 
    }

    for (int n = 1; n <= numOfShells; n++) {//Adds the timestamp to the shells
        workingArray.at(n-1).timestamp = particle.at(0).timestamp;
    }

    // Store working array into the main data array then reset working array
    dataArray.push_back({workingArray});
    workingArray.clear();
    workingArray = std::vector<Shell>(numOfShells, {0,0,0,0});
    return;
}

int main(int argc, char* argv[])  {
    if (argc != 4)  {//Checks to see if the number of arguements is invlaid, produces error message
        std::cout << "Please pass the input and output filename arguments respectively, and then the radius of the shell size for comparison. \nExample: test <input-filename> <out-filename> <shell-size-radius>" << std::endl;
        return 0;      
    }
    
    // Set arguments respectively as defined by user
    std::ifstream InputFile(argv[1]);
    std::ofstream OutputFile(argv[2]);  
    shellWidth = atof(argv[3]);

    // Check if array is empty, initialize if it isn't
    if (workingArray.size() == 0) {
        workingArray = std::vector<Shell>(1, {0,0,0,0});
    }

    while(std::getline(InputFile, line)){
        stream.seekg(0, std::stringstream::beg);
        stream.str(line);
        stream >> inputString;

        if(inputString == "Generated"){//This checks to see if the next line is starting a new timestamp
            stream.seekg(31, std::stringstream::beg);//Skips to where the timestamp will be in the line
            stream >> timestamp;
   
            std::getline(InputFile, line);
            stream.str(line);
            stream.seekg(0, std::stringstream::beg);
            stream >> partNum;

            molNum = (partNum/4);
            if (particle.size() == 0)
                particle = std::vector<Particle>(molNum, {0,0,0,0,0,0,0,0});

            for (int i = 1; i < molNum +1; i++) {

                std::getline(InputFile, line); // Get Info on the Oxygen
                stream.str(line);

                stream.seekg(20, std::stringstream::beg);
                stream >> x >> y>> z;
                xCoM = xCoM + x; // Builds the center of mass vector for the configuration
                yCoM = yCoM + y;
                zCoM = zCoM + z;

                std::getline(InputFile, line); // Get Info on the 1st Hydrogen
                stream.str(line);

                stream.seekg(20, std::stringstream::beg);
                stream >> H1x >> H1y >> H1z;

                std::getline(InputFile, line); // Get Info on the 2nd Hydrogen
                stream.str(line);

                stream.seekg(20, std::stringstream::beg);
                stream >> H2x >> H2y >> H2z;

                std::getline(InputFile, line); // Get Info on the virtual site
                stream.str(line);

                stream.seekg(20, std::stringstream::beg);
                stream >> Vx >> Vy >> Vz;

                xDip = (CHARGE)*(H1x + H2x- 2*Vx);// Find Net Dipole moment of Molecule
                yDip = (CHARGE)*(H1y + H2y- 2*Vy);
                zDip = (CHARGE)*(H1z + H2z- 2*Vz);


                particle.at(i-1) = {timestamp, x, y, z, distCoM, xDip, yDip, zDip}; // Input all relevant info into a particle structure
            }   
        } 
        else {
            xCoM = xCoM/molNum, yCoM = yCoM/molNum, zCoM = zCoM/molNum;
            std::cout << "First pass for timestamp t=" << timestamp << " has completed" <<std::endl;
            CoMAdjustment(OutputFile);
            xCoM = yCoM = zCoM = 0;//Reset CoM variables for the system
        }
    }
    CorrelationCalulation(OutputFile);//Calaculates Correlation function data.
    std::cout << "Done" << std::endl;
return 0;
}
