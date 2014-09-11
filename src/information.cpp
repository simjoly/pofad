//
//  information.cpp
//  
//
//  Created by Simon Joly on 13-11-27.
//  Copyright (c) 2013 __MyCompanyName__. All rights reserved.
//

#include <iostream>
#include <string>
#include <vector>
#include "information.h"

using namespace std;

information::information() {
    this->Verbose = false;
    this->OrganismsFile = false;
    this->DatasetFile = false;
    this->OutputFile = true;
    this->BatchFile = false;
    this->SequenceFile = false;
    this->JukesCantor = true;
    this->IgnoreMissingData = false;
    this->EstimateMissingDist = true;
    this->DistanceMethod = 0;
    this->StandDistanceFileOutput = true;
    this->GapHandling = 0;
    this->output_file = "pofad_output.nex";
}
information::~information(){
	delete [] AllelicInfo;
}
void information::SetVerbose(){
    this->Verbose = true;
}
bool information::IsVerbose(){
    return Verbose;
}
void information::SetIsBatchFile(){
    this->BatchFile = true;
}
bool information::ReturnIsBatchFile(){
    return BatchFile;
}
void information::IsOrganismsFile(){
    this->OrganismsFile = true;
}
bool information::ReturnIsOrganismsFile(){
    return OrganismsFile;
}
void information::IsDatasetFile(){
    this->DatasetFile = true;
}
bool information::ReturnIsDatasetFile(){
    return DatasetFile;
}
void information::IsOutputFile(){
    this->OutputFile = true;
}
bool information::ReturnIsOutputFile(){
    return OutputFile;
}
void information::IsSequenceFile(){
    this->SequenceFile = true;
}
bool information::ReturnIsSequenceFile(){
    return SequenceFile;
}
void information::SetIsJukesCantor(bool abool){
    this->JukesCantor = abool;
}
bool information::ReturnIsJukesCantor(){
    return JukesCantor;
}
void information::SetDistanceMethod(int method) {
    this->DistanceMethod = method;
}	
int information::ReturnDistanceMethod() {
	return DistanceMethod;
}
void information::SetOrganismsFile(string a_file){
    this->file_of_organisms = a_file;
}
string information::ReturnOrganismsFile(){
    return file_of_organisms;
}
void information::SetDatasetsFile(string a_file){
    this->file_of_datasets = a_file;
}
string information::ReturnDatasetsFile(){
    return file_of_datasets;
}
void information::SetOutputFile(string a_file){
    this->output_file = a_file;
}
string information::ReturnOutputFile(){
    return output_file;
}
void information::SetSequenceFile(string a_file){
    this->file_for_sequences = a_file;
}
string information::ReturnSequenceFile(){
    return file_for_sequences;
}
void information::SetGapHandling(int a_number) {
    this->GapHandling = a_number;
}
int information::ReturnGapHandling() {
    return GapHandling;
}
void information::GetDatasetName(string a_string) {
    datasets.push_back(a_string);
}
string information::ReturnDatasetName(int a_number) {
    return datasets[a_number];
}
void information::GetDatasetsNumber(int a_number) {
    this->Number_of_datasets = a_number;
}
int information::ReturnDatasetsNumber() {
    return Number_of_datasets;
}
int information::ReturnDatasetsSize() {
    return datasets.size();
}
void information::SetStandDistanceFileOutput() {
    this->StandDistanceFileOutput = true;
}
void information::SetRawDistanceFileOutput() {
    this->StandDistanceFileOutput = false;
}
bool information::ReturnIsStandDistanceFileOutput() {
    return StandDistanceFileOutput;
}
void information::InitAllelicInfo(int NOrg, int NDatasets) {
    int x,y;
    AllelicInfo = new int* [NOrg];   //Allocate place for matrix
    for (y=0; y < NOrg; y++)
        AllelicInfo[y] = new int[NDatasets];
    if (!AllelicInfo)
    {
        cerr << "Can't allocate space for creating the distance matrix of organisms" << endl;
    }
    for (x=0; x < NOrg; x++) //Initialize the matrix with 0s
    {
        for (y=0; y < NDatasets; y++)
        {
            AllelicInfo[x][y] = 0;
        }
    }
}
void information::EnterNbAllelesforOrganism(int Organism, int dataset, int NbAlleles) {
    AllelicInfo[Organism][dataset] = NbAlleles;
}
int information::ReturnNbAllelesforOrganism(int Organism, int dataset) {
    return AllelicInfo[Organism][dataset];
}

void information::SetIsIgnoreMissingData(bool boolvalue) {
	this->IgnoreMissingData = boolvalue;
}
bool information::ReturnIsIgnoreMissingData() {
	return IgnoreMissingData;
}
void information::SetIsEstimateMissingDist(bool boolvalue) {
	this->EstimateMissingDist = boolvalue;
}
bool information::ReturnIsEstimateMissingDist() {
	return EstimateMissingDist;
}
