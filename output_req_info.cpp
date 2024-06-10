#include <iostream>
#include <fstream>
#include <string>

int main(int argc, char* argv[]) {
    // Check if the filename is provided as a command-line argument
    if (argc != 2) {
        std::cerr << "Usage: " << argv[0] << " <filename>" << std::endl;
        return 1;
    }

    // Open the file
    std::ifstream file(argv[1]);
    if (!file.is_open()) {
        std::cerr << "Error: Unable to open file." << std::endl;
        return 1;
    }

    // Read lines from the file
    std::string line;
    while (std::getline(file, line)) {
        // Check if the line starts with "Total Time taken"
        if (line.find("+ ./kBFS-Ecc -out result.txt inputs/") == 0) {
            std::cout << "Line: " << line << std::endl;
            // Exit the loop after finding the line
        }
        if (line.find("total time excluding") == 0) {
            std::cout << "Line: " << line << std::endl;
            // Exit the loop after finding the line
        }
    }

    // Close the file
    file.close();

    return 0;
}
