.PHONY: setup dev-install clean test

# Setup development environment
setup:
	@echo "🚀 Setting up RunKMC development environment..."
	python3 -m venv .venv || true
	.venv/bin/pip install --upgrade pip
	.venv/bin/pip install build wheel cmake
	.venv/bin/pip install -e .
	python3 build_binary.py
	@echo "✅ Setup complete! Activate with: source .venv/bin/activate"

# Install in development mode
dev-install:
	pip install -e .

# Clean build artifacts
clean:
	rm -rf build/
	rm -rf cpp/build/
	rm -rf dist/
	rm -rf *.egg-info/
	find . -name "__pycache__" -exec rm -rf {} + || true
	find . -name "*.pyc" -delete || true