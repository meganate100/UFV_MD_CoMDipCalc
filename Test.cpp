#include <iostream>
#include <string>
#include <vector>
#include <fstream>
#include <sstream>
#include <cmath>

// Structs
struct Particle {
    float timestamp, x, y, z, distCoM, xDip, yDip, zDip; //This will be where we store the final data for 1 particle. Right now I'm just using it as a general vector.
};

struct Shell {
    float timestamp, NetXDip, NetYDip, NetZDip; //This will be where we store the timestamp and net dipol moment for a given shell.
};

// Variable Declaration

const float CHARGE = 0.6732053036584803;

std::string line;
std::string inputString;
std::istringstream stream(line);
std::vector<Particle> particle;
std::vector<std::vector<Shell>> dataArray;
std::vector<Shell> workingArray;

int partNum;
int molNum;
int currentShellNum;
int numOfShells;
int timeIndex = -1;

float shellWidth;
float xCoM = 0, yCoM = 0, zCoM = 0;
float x, y, z;
float distCoM = 0;
float H1xCoM, H1yCoM, H1zCoM;
float H2xCoM, H2yCoM, H2zCoM;
float VxCoM, VyCoM, VzCoM;
float xDip, yDip, zDip;
float timestamp;
float dt, t0, tf = 0;

//Function that will take a string which describes the particle type (eg: "O" for oxygen, "H" for Hydrogen) and translate it to the appropriate mass
float DotProduct(float x1, float y1, float z1, float x2, float y2, float z2)  {   
    float result = ((x1*x2) + (y1*y2) + (z1*z2));
    return result;
}

// Get max shell size and adjust array size if new shells are required
// Pass center of mass of particle as argument
void CheckShellSize(float CoM) {
    currentShellNum = std::ceil(CoM/shellWidth);
    if (currentShellNum > numOfShells) {
        numOfShells = currentShellNum;
        for (int n = workingArray.size(); n < numOfShells; n++) {
            workingArray.push_back({0,0,0,0});
        }
    }
}

void CorrelationCalulation(std::ofstream& out) {
    t0 = dataArray.at(0).at(0).timestamp;
    dt = dataArray.at(1).at(0).timestamp - t0; //Finds the difference between timesteps
    tf = dataArray.back().at(0).timestamp;
    int n = dataArray.size();
    int k = 0;
    float cor = 0;
    out << "Timesteps: " << dt << "\nInitial Time: " << t0 << "\nFinal Time: " << tf << "\n";
    out << "Number of Shells: " << numOfShells << "\nShell Width: " << shellWidth << "\n";
    for(int i = 0; i < dataArray.back().size() - 1; i++) {
        out << "Correlation for shell: " << i << "\n";
        while (n-k > 0) {
            for(int j = 0; j < n-k; j++){
                cor = cor + (DotProduct(dataArray.at(j).at(i).NetXDip, dataArray.at(j).at(i).NetYDip, dataArray.at(j).at(i).NetZDip,
                                        dataArray.at(j+k).at(i).NetXDip, dataArray.at(j+k).at(i).NetYDip, dataArray.at(j+k).at(i).NetZDip))/(n-k);
            }
            out << "k = " <<  k <<  " n = " <<  n <<"\nCorrelation at t = " << t0 + k*dt << " is " << cor << "\n";
            cor = 0;
            k++;
        }
        k = 0;
    }
}

void CoMAdjustment(std::ofstream& test) {
    //std::cout << "CoM function has been entered" << std::endl;
    for (int i = 1; i < molNum +1; i++) {
        x = particle.at(i-1).x - xCoM;//Get the position vectors of the particles to be relative to the CoM of the system.
        y = particle.at(i-1).y - yCoM;
        z = particle.at(i-1).z - zCoM;

        distCoM = sqrt(DotProduct(x, y, z, x, y, z));//Find shellWidth distance from the CoM of the system.

        CheckShellSize(distCoM);

        particle.at(i-1).x = x;//Re-enters the position values relative to the CoM of the system for all the particles
        particle.at(i-1).y = y;                
        particle.at(i-1).z = z;
        particle.at(i-1).distCoM = distCoM;

        for(int n = 1; n <= numOfShells; n++) {//This creates the net dipole moment of a shellWidth shell of a given timestamp.
            if((distCoM/shellWidth) <= n){
                //std::cout << "The current Net X Dipole moment of Shell #"<< n <<" is: " << workingArray.at(n-1).NetXDip <<
                //" and the next particle's x Dipole moment is: " << particle.at(((timeIndex)*molNum)+i-1).xDip << std::endl;

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
    if (argc != 4)  {
        std::cout << "Please pass the input and output filename arguments respectively, and then the radius of the shell size for comparison. \nExample: test <input-filename> <out-filename> <shell-size-radius>" << std::endl;
        return 0;      
    }
    
    // Set arguments respectively as defined by user
    std::ifstream file(argv[1]);
    std::ofstream TestOutput(argv[2]);  
    shellWidth = atof(argv[3]);

    // Check if array is empty, initialize if it isn't
    if (workingArray.size() == 0) {
        workingArray = std::vector<Shell>(1, {0,0,0,0});
    }

    while(std::getline(file, line)){
        stream.seekg(0, std::stringstream::beg);
        stream.str(line);
        stream >> inputString;

        if(inputString == "Generated"){//This checks to see if the next line is starting a new timestamp.
            timeIndex++;
            stream.seekg(31, std::stringstream::beg);//Skips to where the timestamp will be in the line
            stream >> timestamp;
            //TestOutput << std::endl << timestamp << std::endl;     
            std::getline(file, line);
            stream.str(line);
            stream.seekg(0, std::stringstream::beg);
            stream >> partNum;
            //TestOutput << partNum << std::endl;
            molNum = (partNum/4);
            if (particle.size() == 0)
                particle = std::vector<Particle>(molNum, {0,0,0,0,0,0,0,0});

            for (int i = 1; i < molNum +1; i++) {
                //TestOutput << std::endl << i;
                std::getline(file, line); // Get Info on the Oxygen
                stream.str(line);
                //TestOutput << line << std::endl;
                stream.seekg(20, std::stringstream::beg);
                stream >> x >> y>> z;
                xCoM = xCoM + x; // Builds the center of mass vector for the configuration
                yCoM = yCoM + y;
                zCoM = zCoM + z;

                std::getline(file, line); // Get Info on the 1st Hydrogen
                stream.str(line);
                //TestOutput << line << std::endl;
                stream.seekg(20, std::stringstream::beg);
                stream >> H1xCoM >> H1yCoM >> H1zCoM;

                std::getline(file, line); // Get Info on the 2nd Hydrogen
                stream.str(line);
                //TestOutput << line << std::endl;
                stream.seekg(20, std::stringstream::beg);
                stream >> H2xCoM >> H2yCoM >> H2zCoM;

                std::getline(file, line); // Get Info on the virtual site
                stream.str(line);
                //TestOutput << line << std::endl;
                stream.seekg(20, std::stringstream::beg);
                stream >> VxCoM >> VyCoM >> VzCoM;

                xDip = (CHARGE)*(H1xCoM + H2xCoM- 2*VxCoM);// Find Net Dipole moment of Molecule
                yDip = (CHARGE)*(H1yCoM + H2yCoM- 2*VyCoM);
                zDip = (CHARGE)*(H1zCoM + H2zCoM- 2*VzCoM);


                particle.at(i-1) = {timestamp, x, y, z, distCoM, xDip, yDip, zDip}; // Input all relevant info into a particle structure
                // TestOutput << "Timestamp: " << particle.at(i-1).timestamp << "  x: " << particle.at(i-1).x 
                // << "  y: " << particle.at(i-1).y << "  z: " << particle.at(i-1).z << "  Net xDipole: " 
                // << particle.at(i-1).xDip << "  Net yDipole: " << particle.at(i-1).yDip 
                // << "  Net zDipole: " << particle.at(i-1).zDip << std::endl;
            }   
        } 
        else {
            xCoM = xCoM/molNum, yCoM = yCoM/molNum, zCoM = zCoM/molNum;
            //std::cout << "First pass for timestamp t=" << timestamp << " has completed" <<std::endl;
            //std::cout << "The center of mass for the configuration at t=" << timestamp << " is xCoM: " << xCoM << " yCoM: " << yCoM << " zCoM: " << zCoM << std::endl;
            CoMAdjustment(TestOutput);
            xCoM = yCoM = zCoM = 0;//Reset CoM variables for the system
        }
    }
    CorrelationCalulation(TestOutput);
    std::cout << "Done" << std::endl;
return 0;
}
