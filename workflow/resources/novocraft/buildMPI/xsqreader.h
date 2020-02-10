//
// File:   xsqreader.h
// Author: Colin Hercus
// (C) 2014 Novocraft Technologies Sdn Bhd
//
// Created on 17th March 2014
//

#ifndef _XSQREADER_H
#define	_XSQREADER_H
//#include "ksindex.h"
#include "H5Cpp.h"
#include <string>
using namespace std;


class xsqReader {
    private:
        string xsqFilename;
        int nTags;
        H5::H5File xFile;
        string tagName[2];
        char firstBase[2];
        string tgtLibraryName;
        string crntLibraryName;
        int totalFragments;
        int biggestTile;
        unsigned int readLength[2];
        long iLibrary;                  /// Object index of current library
        H5::Group crntLibrary;          /// Current Library Group
        long iTile;                     /// Object index of current tile
        string tileName;                /// TileName of current tile
        string clippedTileName;
        long nTiles;                    /// Number of Tiles in Library
        long iFragment;                 /// Index of current fragment
        unsigned int nFragments;        /// Number of fragments in the current tile
        unsigned char* colourCallQV[2];       /// Colour Call Data for two tags;
        unsigned short *yxLocation;
        void getTagNames();
        void scanLibraries();
        bool useLibrary(const string &libraryName) const;
        long getIntAttribute(const H5::H5Object &dataset, const char *aName) const;
        void getStringAttribute(const H5::H5Object &ds, const char *aName, string &value) const;

    public:
        xsqReader() 
        {
            colourCallQV[0] = NULL;
            colourCallQV[1] = NULL;
            yxLocation = NULL;
        };
        
        ~xsqReader() 
        {
            xFile.close();
            delete[] colourCallQV[0];
            delete[] colourCallQV[1];
            delete[] yxLocation;
        };
        
        bool open(string &xsqF, string &xsqLibrary);
        bool isPaired() { return nTags > 1; };
        bool next();
        bool eof() { return iLibrary > (long)xFile.getNumObjs(); }
        void read(string &hdr, string &colours, string &qualities, int side);
        unsigned int bp(int side) { return readLength[side]; }
};
#endif  // _XSQREADER_H
