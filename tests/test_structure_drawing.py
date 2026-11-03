"""Structure drawings follow the app theme and sit on a transparent background."""
import io

import pytest

pytest.importorskip("rdkit")
from PIL import Image
from rdkit import Chem

from src.metabolomicsresults import draw_structure

CAFFEINE = "CN1C=NC2=C1C(=O)N(C(=O)N2C)C"
IBUPROFEN = "CC(C)CC1=CC=C(C=C1)C(C)C(=O)O"  # mostly carbon bonds


def _bond_colour(png):
    """Most frequent opaque colour, i.e. the carbon bonds."""
    im = Image.open(io.BytesIO(png)).convert("RGBA")
    opaque = [p[:3] for p in im.getdata() if p[3] == 255]
    return max(set(opaque), key=opaque.count), im


@pytest.mark.parametrize("theme", ["light", "dark"])
def test_background_is_transparent(theme):
    _, im = _bond_colour(draw_structure(Chem.MolFromSmiles(CAFFEINE), theme))
    assert im.getpixel((0, 0))[3] == 0


def test_bonds_contrast_with_theme():
    light, _ = _bond_colour(draw_structure(Chem.MolFromSmiles(IBUPROFEN), "light"))
    dark, _ = _bond_colour(draw_structure(Chem.MolFromSmiles(IBUPROFEN), "dark"))
    assert sum(light) < 3 * 100  # dark bonds on a light page
    assert sum(dark) > 3 * 155  # light bonds on a dark page


def test_unparsable_structure_draws_nothing():
    assert draw_structure(Chem.MolFromSmiles("not a smiles")) is None
