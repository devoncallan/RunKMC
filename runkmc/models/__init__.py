from typing import Dict, Any
from pathlib import Path


from runkmc import PATHS

REGISTERED_MODELS = {
    "FRP2": PATHS.TEMPLATE_DIR / "binary/FRP2_Template.txt",
    "CRP1": PATHS.TEMPLATE_DIR / "binary/CRP1_Template.txt",
    "CRP3": PATHS.TEMPLATE_DIR / "binary/CRP3_Template.txt",
}


def create_input_file(
    model_name: str, kmc_inputs: Dict[str, Any], filepath: Path | str
) -> None:

    model_filepath = REGISTERED_MODELS.get(model_name)
    if model_filepath is None:
        raise ValueError(f"Model {model_name} not supported")

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
