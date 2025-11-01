from __future__ import annotations
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Any, Optional

TEMPLATES_DIR = Path(__file__).parent

REGISTERED_TEMPLATES = {
    "FRP1": TEMPLATES_DIR / "homo/FRP1.template",
    "BUP2": TEMPLATES_DIR / "binary/BUP2.template",
    "BUP2_Rev": TEMPLATES_DIR / "binary/BUP2_Rev.template",
    "FRP2": TEMPLATES_DIR / "binary/FRP2.template",
    "CRP1": TEMPLATES_DIR / "binary/CRP1.template",
    "CRP3": TEMPLATES_DIR / "binary/CRP3.template",
}


@dataclass
class TemplateParameter:
    name: str
    default: Optional[str] = None
    required: bool = True

    def __post_init__(self):
        self.required = self.default is None


class Template:
    """Input file template with parameter validation."""

    PARAM_PATTERN = re.compile(r"\{\{(\w+)(?:=([^}]+))?\}\}")

    def __init__(self, content: str):
        self.content = content
        self.parameters = self._parse_parameters()

    def _parse_parameters(self) -> Dict[str, TemplateParameter]:
        params = {}
        for match in self.PARAM_PATTERN.finditer(self.content):
            name = match.group(1)
            default = match.group(2)
            if name not in params:
                params[name] = TemplateParameter(name, default)
        return params

    @property
    def required_params(self) -> list[str]:
        return [p.name for p in self.parameters.values() if p.required]

    @property
    def optional_params(self) -> Dict[str, str]:
        return {
            p.name: p.default
            for p in self.parameters.values()
            if not p.required and p.default
        }

    def get_info(self) -> Dict[str, Any]:
        return {
            "required": self.required_params,
            "optional": self.optional_params,
            "all": list(self.parameters.keys()),
        }

    def validate(self, user_params: Dict[str, Any]) -> Dict[str, str]:
        missing = set(self.required_params) - set(user_params.keys())
        if missing:
            raise ValueError(f"Missing required parameters: {sorted(missing)}")

        unknown = set(user_params.keys()) - set(self.parameters.keys())
        if unknown:
            raise ValueError(f"Unknown parameters: {sorted(unknown)}")

        complete_params = {}
        for name, param in self.parameters.items():
            if name in user_params:
                complete_params[name] = str(user_params[name])
            elif param.default is not None:
                print(
                    f"\tParameter `{name}` not provided. Using default: {param.default}"
                )
                complete_params[name] = param.default
            else:
                raise ValueError(f"Parameter {name} has no value and no default")

        return complete_params

    def render(self, user_params: Dict[str, Any]) -> str:
        complete_params = self.validate(user_params)

        def replacer(match):
            param_name = match.group(1)
            return complete_params[param_name]

        return self.PARAM_PATTERN.sub(replacer, self.content)

    @staticmethod
    def from_file(filepath: Path | str) -> Template:
        filepath = Path(filepath)
        if not filepath.exists():
            raise FileNotFoundError(f"Template not found: {str(filepath)}")
        with open(filepath, "r") as f:
            return Template(f.read())

    @staticmethod
    def load_registered(template_name: str) -> Template:
        if isinstance(template_name, str) and template_name in REGISTERED_TEMPLATES:
            template_path = REGISTERED_TEMPLATES[template_name]
        else:
            raise ValueError(f"Template not registered: {template_name}")

        return Template.from_file(template_path)

    @staticmethod
    def load(filepath_or_name: Path | str) -> Template:
        try:
            return Template.load_registered(str(filepath_or_name))
        except ValueError:
            return Template.from_file(filepath_or_name)
        except Exception as e:
            raise ValueError(f"Could not not find or load template: {filepath_or_name}")


def create_input_file(
    template_name: str | Path,
    parameters: Dict[str, Any],
    filepath: Path | str,
    validate: bool = True,
) -> None:

    template = Template.load(template_name)

    if validate:
        content = template.render(parameters)
    else:
        content = template.content
        for key, value in parameters.items():
            content = content.replace(f"{{{{{key}}}}}", str(value))

    filepath = Path(filepath)
    filepath.parent.mkdir(parents=True, exist_ok=True)
    with open(filepath, "w") as f:
        f.write(content)
