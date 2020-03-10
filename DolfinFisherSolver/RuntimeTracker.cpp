
#include "RuntimeTracker.h"

// class RuntimeInformation
///////////////////////////////
RuntimeTracker::RuntimeTracker(){}
RuntimeTracker::RuntimeTracker(int verbose, bool toCsv, std::string csvPath){
	// general state
	verbose_ = verbose;					// verbosity of tracker
	elapsedAll_ = -1;					// if elapsedAll_ > 0, tracking has been ended
	// iterations state
	inIteration_ = false;				// no running iteration currently
	numberIterations_ = 0;				// number of tracked iterations
	// output of iteration data
	toCsv_ = toCsv;						// if tracker outputs to given csv file
	csvPath_ = csvPath;					// path to file storing iteration data as cvs
	if(toCsv_){
		csv_.open(csvPath, std::ios_base::trunc);
		csv_ << "iteration data file \n";
		csv_.close();
	}
}
std::string RuntimeTracker::asString(){
	std::stringstream ss;
	// tracking has been ended, add simulation recap
	if(elapsedAll_ > 0){
		int minutes = elapsedAll_ / 60;
		ss << std::endl << "Simulation finished:\nnumber iterations: "<< numberIterations_ << std::endl <<
				"elapsed time: " << (minutes / 60) << "h " << (minutes % 60) << " min" << std::endl;
	}
	return ss.str();
}
// start/end of a simulation
void RuntimeTracker::startTracking(){
	startTimeRun_ = std::chrono::system_clock::now();
	// print
	if(verbose_ > 1) {
		std::time_t now = std::chrono::system_clock::to_time_t(startTimeRun_);
		std::string time(30, '\0');
		std::strftime(&time[0], time.size(), "%Y-%m-%d %H:%M:%S", std::localtime(&now));
		std::cout << "start tracker: (" << time << ")" << std::endl;
	}
}
void RuntimeTracker::stopTracking(){
	endTimeRun_ = std::chrono::system_clock::now();
	std::chrono::duration<double> duration = endTimeRun_ - startTimeRun_;
	elapsedAll_ = duration.count();
	// print
	if(verbose_ > 1) { // verbose level 2
		std::time_t now = std::chrono::system_clock::to_time_t(endTimeRun_);
		std::string time(30, '\0');
		std::strftime(&time[0], time.size(), "%Y-%m-%d %H:%M:%S", std::localtime(&now));
		std::cout << "stop tracker: (" << time << ")" << std::endl <<
				"	elapsed time: " << elapsedAll_/60 << " min" << std::endl <<
				"	elapsed time: " << elapsedAll_/3600 << " h" << std::endl;
	}
}
//  start/end of a new iteration
void RuntimeTracker::newIteration(){
	inIteration_ = true;
	std::vector<double> empty = {};
	currIteration_ = empty;
}
void RuntimeTracker::startTime(){
	startTimeIteration_ = std::chrono::system_clock::now();
}
void RuntimeTracker::endTime(){
	auto endTimeIteration = std::chrono::system_clock::now();
	std::chrono::duration<double> duration = endTimeIteration - startTimeIteration_;
	elapsedIt_ = duration.count();
}
void RuntimeTracker::endIteration(){
	if(inIteration_){
		inIteration_ = false;
		currIteration_.push_back(elapsedIt_);
		iterations_.push_back(currIteration_);
		numberIterations_++;
		// print
		if(verbose_ > 2){ // verbose level 3
			std::cout << "iteration " << numberIterations_ << " at t= " << currIteration_.at(0) << std::endl <<
					"	dt = " << currIteration_.at(1) << std::endl <<
					"	residual = " << currIteration_.at(2) << std::endl;
			std::cout << std::endl;
		}
	}
	else{
		std::cout << "WARNING: seems like no iteration has been started yet!" << std::endl <<
				"	No information stored" << std::endl;
	}

	// if buffer limit is reached: write iteration data to csv file and clear buffer
	if(toCsv_ && numberIterations_ % iterationBufferSize_ == 0){
		if(verbose_ > 2){	// verbosity level 3
			std::cout << "writing iteration data (last " << iterationBufferSize_ << " of " << numberIterations_ << ") to csv..." << std::endl << std::endl;
		}

		int numIts = iterations_.size();
		csv_.open(csvPath_, std::ios_base::app);
		for(int i = 0; i < numIts; i++){
			auto it = iterations_[i];
			for(int j = 0; j < it.size(); j++){
				csv_ << it.at(j) << ",";
			}
			csv_ << "\n";
		}
		csv_.close();
		iterations_.clear();
	}
}
void RuntimeTracker::addIterationFormat(std::string format){
	format_ = format;
	if(toCsv_){
		csv_.open(csvPath_, std::ios_base::trunc);
		csv_ << format << ", elapsed time \n";
		csv_.close();
	}
}
// add data from newton solver to current iteration instance
void RuntimeTracker::addIterationData(std::vector<double> data){
	currIteration_ = data;
}


