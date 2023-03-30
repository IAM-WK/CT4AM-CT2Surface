#pragma once
#include <chrono>

namespace progress {
	// for progressbar
	// ready to called in parallel with multiple instances. Caller has to make sure this is called in OMP critical sections
	static float processed = -1, oldprocessed = -1; // for progress bar
	auto TimeBefore = std::chrono::steady_clock::now();
	auto TimeAfter = std::chrono::steady_clock::now();
	auto DeltaSum = std::chrono::duration_cast<std::chrono::milliseconds>(TimeAfter - TimeBefore).count();
	
	// to be called for displaying progress in percent
	void progressbar(const size_t& progresscounter, const size_t& ToBeDone) {
		oldprocessed = processed;
		processed = std::roundf(static_cast<float>(progresscounter) / static_cast<float>(ToBeDone) * 100); // processed in percent (0.xy)
		if (oldprocessed != processed) { // new percent values only
			std::cout << "\rprocessed " << processed << "% of projections!";
		}
	}

	// to be called for displaying progress in percent and time
	// ofc only works in serial execution...
	void progressbarTime(const size_t& progresscounter, const size_t& ToBeDone) {
		oldprocessed = processed;
		processed = std::roundf(static_cast<float>(progresscounter) / static_cast<float>(ToBeDone) * 100); // processed in percent (0.xy)

		if (oldprocessed != processed) { // new percent values only
			TimeAfter = std::chrono::steady_clock::now();
			DeltaSum += std::chrono::duration_cast<std::chrono::milliseconds>(TimeAfter - TimeBefore).count();
			size_t TimeLeft = static_cast<size_t>(std::roundf(((DeltaSum / processed) * (100.0) - DeltaSum) / 1000.0));
			size_t TimeLeftInMinutes = static_cast<size_t>(std::floor(static_cast<size_t>(TimeLeft) / 60.0));
			bool ShowinMinutes = TimeLeft > 60.0;

			std::cout << "\rprocessed " << processed << "% of projections! ";
			if (progresscounter == 0) {
				return;
			}
			//std::cout << "Estimated time left: " << std::to_string((std::chrono::duration_cast<std::chrono::milliseconds>(TimeAfter - TimeBefore).count()) * (100.0) / 1000.0) << " seconds";
			if (ShowinMinutes) {
				std::cout << "Estimated time left: " << std::to_string(TimeLeftInMinutes) << "min" << std::to_string(static_cast<size_t>(TimeLeft - (TimeLeftInMinutes * 60.0))) << "sec                 ";
			}
			else {
				std::cout << "Estimated time left: " << std::to_string(static_cast<size_t>(std::roundf(((DeltaSum / processed) * (100.0) - DeltaSum) / 1000.0))) << " seconds                      ";
			}
			TimeBefore = TimeAfter;
		}
	}
}
