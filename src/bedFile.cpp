/*****************************************************************************
bedFile.cpp

(c) 2009 - Aaron Quinlan
Hall Laboratory
Department of Biochemistry and Molecular Genetics
University of Virginia
aaronquinlan@gmail.com

Licensed under the GNU General Public License 2.0 license.
******************************************************************************/
#include "bedFile.h"

/*******************************************
Class methods
*******************************************/

// Constructor
BedFile::BedFile(string bedFile)
: bedFile(bedFile),
_typeIsKnown(false),
_lineNum(0)
{}

// Destructor
BedFile::~BedFile(void) {
}


void BedFile::Open(void) {
    if (bedFile == "stdin") {
        _bedStream = &cin;
    }
    // New method thanks to Assaf Gordon
    else if ((isGzipFile(bedFile) == false) && (isRegularFile(bedFile) == true)) {
       // open an ifstream
        ifstream beds(bedFile.c_str(), ios::in);

        // can we open the file?
        if ( !beds ) {
            cerr << "Error: The requested bed file (" << bedFile << ") could not be opened. Exiting!" << endl;
            exit (1);
        }
        else {
            // if so, close it (this was just a test)
            beds.close();       
            // now set a pointer to the stream so that we
            _bedStream = new ifstream(bedFile.c_str(), ios::in);
        }
    } 
    else if ((isGzipFile(bedFile) == true) && (isRegularFile(bedFile) == true)) {        
        igzstream beds(bedFile.c_str(), ios::in);
        if ( !beds ) {
            cerr << "Error: The requested bed file (" << bedFile << ") could not be opened. Exiting!" << endl;
            exit (1);
        }
        else {
            // if so, close it (this was just a test)
            beds.close();       
            // now set a pointer to the stream so that we
            _bedStream = new igzstream(bedFile.c_str(), ios::in);
        }
    }
    else {
        cerr << "Error: Unexpected file type (" << bedFile << "). Exiting!" << endl;
        exit(1);
    }
}


// Close the BED file
void BedFile::Close(void) {
    if (bedFile != "stdin") delete _bedStream;
}


BED BedFile::GetNextBed() {

    BED bed;
    
    // make sure there are still lines to process.
    // if so, tokenize, validate and return the BED entry.
    if (_bedStream->good()) {
        string bedLine;
        vector<string> bedFields;
        bedFields.reserve(12);

        // parse the bedStream pointer
        getline(*_bedStream, bedLine);
        _lineNum++;

        // split into a string vector.
        Tokenize(bedLine,bedFields);

        // load the BED struct as long as it's a valid BED entry.
        bed.status = parseLine(bed, bedFields);
        return bed;
    }
    else {
        // default if file is closed or EOF
        bed.status = BED_INVALID;
        return bed;
    }
}

vector<BED> BedFile::FindOverlapsPerBin(string chrom, CHRPOS start, CHRPOS end, float overlapFraction) {
    vector<BED> hits;

    BIN startBin, endBin;
    startBin = (start >> _binFirstShift);
    endBin = ((end-1) >> _binFirstShift);

    // loop through each bin "level" in the binning hierarchy
    for (BINLEVEL i = 0; i < _binLevels; ++i) {

        // loop through each bin at this level of the hierarchy
        BIN offset = _binOffsetsExtended[i];
        for (BIN j = (startBin+offset); j <= (endBin+offset); ++j)  {

            // loop through each feature in this chrom/bin and see if it overlaps
            // with the feature that was passed in.  if so, add the feature to 
            // the list of hits.
            vector<BED>::const_iterator bedItr = bedMap[chrom][j].begin();
            vector<BED>::const_iterator bedEnd = bedMap[chrom][j].end();

            for (; bedItr != bedEnd; ++bedItr) {
                // do we have sufficient overlap?
                float size = (float) end-start;
                if ( (float) overlaps(bedItr->start, bedItr->end, start, end) / size > overlapFraction) {
                        hits.push_back(*bedItr);
                }
            }
        }
        startBin >>= _binNextShift;
        endBin >>= _binNextShift;
    }
    return hits;
}

vector<BED> BedFile::FindOverlapsPerBin(string chrom, CHRPOS start, CHRPOS end, string strand, float overlapFraction) {
    vector<BED> hits;

    BIN startBin, endBin;
    startBin = (start >> _binFirstShift);
    endBin = ((end-1) >> _binFirstShift);

    // loop through each bin "level" in the binning hierarchy
    for (BINLEVEL i = 0; i < _binLevels; ++i) {

        // loop through each bin at this level of the hierarchy
        BIN offset = _binOffsetsExtended[i];
        for (BIN j = (startBin+offset); j <= (endBin+offset); ++j)  {

            // loop through each feature in this chrom/bin and see if it overlaps
            // with the feature that was passed in.  if so, add the feature to 
            // the list of hits.
            vector<BED>::const_iterator bedItr = bedMap[chrom][j].begin();
            vector<BED>::const_iterator bedEnd = bedMap[chrom][j].end();

            for (; bedItr != bedEnd; ++bedItr) {
                // do we have sufficient overlap?
                float size = (float) end-start;
                if ( (float) overlaps(bedItr->start, bedItr->end, start, end) / size > overlapFraction
                    &&
                    (strand == bedItr->strand))
                {
                    hits.push_back(*bedItr);
                }
            }
        }
        startBin >>= _binNextShift;
        endBin >>= _binNextShift;
    }
    return hits;
}


void BedFile::setGff (bool gff) {
    if (gff == true) this->_isGff = true;
    else this->_isGff = false;
}


void BedFile::setVcf (bool vcf) {
    if (vcf == true) this->_isVcf = true;
    else this->_isVcf = false;
}


void BedFile::setFileType (FileType type) {
    _fileType    = type;
    _typeIsKnown = true;
}


void BedFile::setBedType (int colNums) {
    bedType = colNums;
}

void BedFile::loadBedFileIntoMap() {

    BED bed, nullBed;
    BedLineStatus bedStatus;

    Open();
    bed = GetNextBed();
    while ( bed.status != BED_INVALID) {
        if (bed.status == BED_VALID) {
            BIN bin = getBin(bed.start, bed.end);
            bedMap[bed.chrom][bin].push_back(bed);
            bed = nullBed;
        }
        bed = GetNextBed();
    }
    Close();
}
