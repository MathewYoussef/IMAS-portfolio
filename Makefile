.PHONY: quickstart figures audit docs

quickstart:
	@echo "TODO: Execute sample notebooks on lightweight datasets (Agent 3 to supply)."

figures:
	@echo "TODO: Regenerate all figures into scaffold/*/figures once notebooks are wired."

audit:
	@echo "TODO: Run claim checks and health tests across all blocks."

docs:
	python3 -m mkdocs build -q
