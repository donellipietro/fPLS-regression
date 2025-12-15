# Define the Rscript command
RSCRIPT := Rscript

# Targets
.PHONY: install build install_fdaPDE2 install_femR clean distclean tests test_example

# Description for targets

# Default target
all: build

# Installation targets
install_fdaPDE2: 
	@echo "Installing fdaPDE2..."
	@$(RSCRIPT) src/installation/install_fdaPDE2.R
	
install_femR: 
	@echo "Installing femR..."
	@$(RSCRIPT) src/installation/install_femR.R

install: install_fdaPDE2 install_femR
	@echo "Installation completed."

# Test targets

test_1D: build clean
	@echo "Running fPLS - 1D tests..."
	# @./run_tests_parallel.sh 1D test1
	# @./run_tests_parallel.sh 1D test2
	@./run_tests_parallel.sh 1D test3
	
test_2D: build clean
	@echo "Running fPLS - 2D tests..."
	@./run_tests_parallel.sh 2D test1
	# @./run_tests_parallel.sh 2D test2
	# @./run_tests_parallel.sh 2D test3
	
tests: test_2D
	@echo "All tests completed."

# Build target
build:
	@echo "Creating necessary directories..."
	@mkdir -p data
	@mkdir -p results
	@mkdir -p images
	@echo "Build completed."

# Clean targets
clean_options:
	@$(RM) -r queue/
	@$(RM) -r logs/


clean: clean_options
	@echo "Cleaning temporary files..."
	@$(RM) *.aux *.log *.pdf *.txt *.json
	@$(RM) .Rhistory
	@$(RM) .RData
	@echo "Cleanup completed."
	
	
distclean: clean
	@echo "Attention! This will remove additional generated files."
	@read -p "Are you sure you want to continue? [y/n]: " confirm && [ "$$confirm" = "y" ] || (echo "Cleanup aborted." && false)
	@echo "Removing additional generated files..."
	@$(RM) -r images/
	@$(RM) -r results/
	@echo "Additional cleanup completed."
