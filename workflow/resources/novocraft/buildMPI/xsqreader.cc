//
// File:   xsqreader.cc
// Author: Colin Hercus
// (C) 2014 Novocraft Technologies Sdn Bhd
//
// Created on 17th March 2014
//

#include <algorithm>
#include "xsqreader.h"
#include <vector>
#include <string.h>
#include <assert.h>

using namespace std;

static const char defaultLibrary[] = "DefaultLibrary";

bool xsqReader::open(string &xsqF, string &xsqLibrary)
{
    if(!xFile.isHdf5(xsqF.c_str())) {
        fprintf(stderr,"XSQ Error: '%s' is not a HDF5 format file.\n", xsqF.c_str());
        exit(-1);        
    }
    try
    {
        xFile.openFile(xsqF.c_str(), H5F_ACC_RDONLY);
    }
    catch(...)
    {
        fprintf(stderr,"XSQ Error: Exception opening '%s'.\n", xsqF.c_str());
        exit(-1);                
    }
    xsqFilename = xsqF;
    getTagNames();
    tgtLibraryName = xsqLibrary;
//    if(tgtLibraryName != "")
//        tgtLibraryName += "_";
    scanLibraries();
    if(totalFragments == 0)
    {
        if(xsqLibrary == "")
            fprintf(stderr,"XSQ Error: No Fragments found for '%s'.\n", defaultLibrary);
        else
            fprintf(stderr,"XSQ Error: No Fragments found for library '%s'.\n", xsqLibrary.c_str());
        exit(-1);            
    }
    iLibrary = -1;
    nTiles = 0;
    iTile = -1;
    nFragments = 0;
    iFragment = -1;
    colourCallQV[0] = new unsigned char[biggestTile*readLength[0]];
    if(isPaired())
        colourCallQV[1] = new unsigned char[biggestTile*readLength[1]];
    yxLocation = new unsigned short [biggestTile * 2];
    return true;
}

bool xsqReader::useLibrary(const string &libraryName) const
{
    if(tgtLibraryName == "" && libraryName == defaultLibrary)
        return true;
    if(libraryName == tgtLibraryName)
        return true;
    if(libraryName.length() < tgtLibraryName.length())
        return false;
    if(tgtLibraryName.compare(0, tgtLibraryName.length(), libraryName, 0, tgtLibraryName.length()) != 0)
        return false;
//    size_t i = libraryName.find_last_of("_");
    return libraryName.find_last_of("_") == tgtLibraryName.length();
}

long xsqReader::getIntAttribute(const H5::H5Object &ds, const char *aName) const
{
    long i = 0;
    const H5::Attribute& at = ds.openAttribute(aName);
    const H5::DataType& dt = at.getDataType();
    at.read(dt, &i); 
    return i;
}

void xsqReader::getStringAttribute(const H5::H5Object &ds, const char *aName, string &value) const
{
    const H5::Attribute& at = ds.openAttribute(aName);
    const H5::DataType& dt = at.getDataType();
    at.read(dt, value); 
}

void xsqReader::getTagNames() 
{
    try
    { 
	const H5::Group &runMetadata = xFile.openGroup("RunMetadata");
	const H5::Group &tagDetails = runMetadata.openGroup("TagDetails");

	int nTagObj = (int)tagDetails.getNumObjs();
        nTags = 0;
	std::vector<std::string> used_tags_names(nTagObj);
        struct ce {
            char dataSetName[12];
            char encoding[3];
            unsigned char stride;
            unsigned int numColorCalls;
        } *ceData;

	for(unsigned t = 0; t < nTagObj; t++) {
            tagName[nTags] = tagDetails.getObjnameByIdx(t);
            const H5::Group &tagGroup = tagDetails.openGroup(tagName[nTags]);
            if(getIntAttribute(tagGroup,"IsColorPresent") != 1)
                continue;
            const H5::DataSet &colourEncoding = tagGroup.openDataSet("DatasetColorEncoding");
            H5::CompType ceType(sizeof(ce));
            H5::StrType dsnType(H5::PredType::C_S1,12); 
            H5::StrType encodingType(H5::PredType::C_S1,3); 
            ceType.insertMember("DataSetName", HOFFSET(ce, dataSetName), dsnType);
            ceType.insertMember("Encoding", HOFFSET(ce, encoding), encodingType);
            ceType.insertMember("Stride", HOFFSET(ce, stride), H5::PredType::NATIVE_UCHAR);
            ceType.insertMember("NumColorCalls", HOFFSET(ce, numColorCalls), H5::PredType::NATIVE_UINT);
            hsize_t rows[1];
            colourEncoding.getSpace().getSimpleExtentDims(rows, NULL);
            ceData = new ce[rows[0]];
            colourEncoding.read(ceData, ceType);
            if(ceData[0].stride != 1)
                continue;
            if(strcmp(ceData[0].encoding, "11") != 0)
                continue;
            readLength[nTags] = ceData[0].numColorCalls;
            std::string tagSequence;
            getStringAttribute(tagGroup, "TagSequence", tagSequence);
            firstBase[nTags] = tagSequence[tagSequence.length() - 1];
            nTags++;
        }
        
        if(nTags == 0 || nTags > 2)
        {
            fprintf(stderr,"XSQ Error: XSQ files has %d tags with Stride 1 and '11' encoding. I can't handle that!", nTags);
            exit(-1);               
        }
    } 
    catch (...) 
    { 
        fprintf(stderr,"XSQ Error: HDF5 Exception in %s:%s:%d\n", __FILE__, __func__, __LINE__);
        exit(-1); 
    } 
}

bool xsqReader::next()
{
    while(true)
    {
        if(eof())
            return false;
        if(++iFragment < nFragments)
            return true;
        if(++iTile < nTiles)
        {
            const H5::Group& crntLibrary = xFile.openGroup(crntLibraryName);
            clippedTileName = tileName = crntLibrary.getObjnameByIdx(iTile);
            size_t nz = tileName.find_first_not_of("0");
            if(nz == string::npos)
                clippedTileName.resize(1);
            else if(nz > 0)
                clippedTileName.replace(0, nz, "");
            const H5::Group &tile = crntLibrary.openGroup(tileName);
            const H5::Attribute& fragmentCountAttribute = tile.openAttribute("FragmentCount");
            const H5::DataType& fragmentCountType = fragmentCountAttribute.getDataType();
            unsigned int fragmentCount;
            fragmentCountAttribute.read(fragmentCountType, &fragmentCount); 
            nFragments = fragmentCount;      
            iFragment = -1;
            {
                const H5::Group &tag = tile.openGroup(tagName[0]);
                const H5::DataSet &dsColorCallQV = tag.openDataSet("ColorCallQV");
                dsColorCallQV.read(colourCallQV[0], H5::PredType::NATIVE_UINT8);
            }
            if(isPaired())
            {
                const H5::Group tag = tile.openGroup(tagName[1]);
                const H5::DataSet &dsColorCallQV = tag.openDataSet("ColorCallQV");
                dsColorCallQV.read(colourCallQV[1], H5::PredType::NATIVE_UINT8);                
            }
            {
            	const H5::Group& fragments = tile.openGroup("Fragments");
                const H5::DataSet& dsYxLocation = fragments.openDataSet("yxLocation");
                dsYxLocation.read(yxLocation, H5::PredType::NATIVE_UINT16);
            }
            continue;
        }
        while(++iLibrary < xFile.getNumObjs())
        {
            crntLibraryName = xFile.getObjnameByIdx(iLibrary);
            if(!useLibrary(crntLibraryName))
                continue;
            nTiles = xFile.openGroup(crntLibraryName).getNumObjs();
            iTile = -1;
            nFragments = 0;
            iFragment = -1;
            break;
        }
    }    
}

void xsqReader::scanLibraries()
{
    biggestTile = 0;
    totalFragments = 0;
    try
    {
	for(unsigned l = 0; l < xFile.getNumObjs(); l++)
	{
//            const std::string& 
            crntLibraryName = xFile.getObjnameByIdx(l);
            if(useLibrary(crntLibraryName))
            {
                const H5::Group &library = xFile.openGroup(crntLibraryName);
                for(unsigned t = 0; t < library.getNumObjs(); t++)
                {
                    const std::string& tileName = library.getObjnameByIdx(t);
                    const H5::Group &tile = library.openGroup(tileName);
                    const H5::Attribute& fragmentCountAttribute = tile.openAttribute("FragmentCount");
                    const H5::DataType& fragmentCountType = fragmentCountAttribute.getDataType();
                    unsigned int fragmentCount;
                    fragmentCountAttribute.read(fragmentCountType, &fragmentCount);
                    if(fragmentCount > biggestTile)
                        biggestTile = fragmentCount;
                    totalFragments += fragmentCount;
                }    
            }
	}                
    }
    catch (...) 
    { 
        fprintf(stderr,"XSQ Error: HDF5 Exception in %s:%s:%d\n", __FILE__, __func__, __LINE__);
        exit(-1); 
    } 
}

void xsqReader::read(string &hdr, string &colours, string &qualities, int side)
{
    /*
     *  Call = CallQV & 0x03;
     *  QV = CallQV >> 2; // 0, 1, 2, and 63 are special
     *   0 and 1 are reserved for future usage; 
     *   2 indicates a missing quality value (not trimming)
     *   63 refers to an uninformative or missing call.
     * 
     *  hdr  ><tileName>_<Y>_<X>_<tagName>
     */
    colours.resize(readLength[side]+1);
    qualities.resize(readLength[side]+1);
    colours[0] = firstBase[side];
    qualities[0] = '!';
    char yx[14]; //Long enough for unsigned int type
    assert(sprintf(yx, "_%u_%u_", yxLocation[iFragment * 2], yxLocation[iFragment * 2 + 1] ) < 14);
    hdr = '@';
    hdr += clippedTileName;
    hdr += yx;
    hdr += tagName[side];
    
    unsigned char* read = colourCallQV[side] + iFragment * readLength[side];
    int i;
    for(i = 0; i < readLength[side]; i++)
    { 
        int q = read[i] >> 2;
        if(q > 2 & q < 63) 
        {
            qualities[i+1] = q + '!';
            colours[i+1] = (read[i] & 0x3) + '0';
        }
        else
        {
            qualities[i+1] = '!';
            colours[i+1] = '.';            
        }
            
    }  
    qualities[i+1] = '\0';
    colours[i+1] = '\0';
}


