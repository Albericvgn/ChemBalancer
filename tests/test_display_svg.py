# tests/test_display_svg.py

from chembalancer.chembalancer import display_svg
import base64

def test_display_svg():
    # Test case 1: Display SVG with a simple circle
    svg_content_1 = """
    <svg height="100" width="100">
      <circle cx="50" cy="50" r="40" stroke="black" stroke-width="3" fill="red" />
    </svg>
    """
    expected_html_1 = f'<img src=\'data:image/svg+xml;base64,{base64.b64encode(svg_content_1.encode("utf-8")).decode("utf-8")}\'/>'
    assert display_svg(svg_content_1) == expected_html_1
    
    # Test case 2: Display SVG with a square
    svg_content_2 = """
    <svg height="100" width="100">
      <rect width="100" height="100" style="fill:rgb(0,0,255);stroke-width:1;stroke:rgb(0,0,0)" />
    </svg>
    """
    expected_html_2 = f'<img src=\'data:image/svg+xml;base64,{base64.b64encode(svg_content_2.encode("utf-8")).decode("utf-8")}\'/>'
    assert display_svg(svg_content_2) == expected_html_2

    # Test case 3: Display SVG with a line
    svg_content_3 = """
    <svg height="100" width="100">
      <line x1="0" y1="0" x2="100" y2="100" style="stroke:rgb(255,0,0);stroke-width:2" />
    </svg>
    """
    expected_html_3 = f'<img src=\'data:image/svg+xml;base64,{base64.b64encode(svg_content_3.encode("utf-8")).decode("utf-8")}\'/>'
    assert display_svg(svg_content_3) == expected_html_3
