from typing import Dict, Any
from pathlib import Path


from runkmc import PATHS

REGISTERED_MODELS = {
    "FRP1": PATHS.TEMPLATE_DIR / "homo/FRP1_Template.txt",
    "FRP2": PATHS.TEMPLATE_DIR / "binary/FRP2_Template.txt",
    "CRP1": PATHS.TEMPLATE_DIR / "binary/CRP1_Template.txt",
    "CRP3": PATHS.TEMPLATE_DIR / "binary/CRP3_Template.txt",
}


def create_input_file(
    model_name: str | Path, kmc_inputs: Dict[str, Any], filepath: Path | str
) -> None:

    if isinstance(model_name, str) and model_name in REGISTERED_MODELS:
        model_filepath = REGISTERED_MODELS.get(model_name)
    elif Path(model_name).exists():
        model_filepath = Path(model_name)
    else:
        raise ValueError(f"Model {model_name} not supported")
    
    if not model_filepath or not model_filepath.exists():
        raise FileNotFoundError(f"Model template file {model_filepath} not found.")

    # Read model template file
    with open(model_filepath, "r") as file:
        template_content = file.read()

    # Replace placeholders in template with input values
    for key, value in kmc_inputs.items():
        placeholder = "{" + key + "}"

        try:
            value = str(value)
        except ValueError:
            raise ValueError(f"Value for {key} cannot be converted to string.")

        template_content = template_content.replace(placeholder, value)

    with open(filepath, "w") as file:
        file.write(template_content)
