#!/usr/bin/env python3
# type: ignore

import sys
import argparse
from typing import Iterator, Any
try:
    import fitz  # PyMuPDF
except ImportError:
    print("The combination of PDF pages for the Feynman Diagram Drawing requires the python module fitz. Please install it with 'python -m pip install PyMuPDF'.")
    sys.exit(1)


def divide_chunks(l: list[Any], n: int) -> Iterator[Any]:
    for i in range(0, len(l), n):
        yield l[i:i + n]


def merge_pages_into_grid(input_pdf_path: str, output_pdf_path: str, pages_to_merge: list[int] | None, grid_shape: tuple[int, int], margins: tuple[float]):
    rows, cols = grid_shape
    pdf = fitz.open(input_pdf_path)

    if pages_to_merge is None:
        pages_to_merge = list(range(len(list(pdf.pages()))))

    # Create a new PDF
    pdf_output = fitz.open()

    # Get dimensions from the first page
    page: list[Any] = pdf[pages_to_merge[0]]
    rect = page.rect
    new_width = rect.width * cols
    new_height = rect.height * rows

    for pages_chunk_to_merge in divide_chunks(pages_to_merge, rows * cols):
        # Create a blank page with new dimensions
        new_page = pdf_output.new_page(width=new_width, height=new_height)

        for idx, page_number in enumerate(pages_chunk_to_merge):
            row = idx // cols
            col = idx % cols
            cell_height = rect.height * (1 - margins[1])
            cell_width = rect.width * (1 - margins[0])
            x_offset = col * rect.width + rect.width * margins[0] / 2
            y_offset = row * rect.height + rect.height * margins[1] / 2
            target_rect = fitz.Rect(
                x_offset, y_offset, x_offset + cell_width, y_offset + cell_height)

            # Insert the page into the new page
            new_page.show_pdf_page(target_rect, pdf, page_number)

    # Save the new PDF
    pdf_output.save(output_pdf_path)


if __name__ == "__main__":

    parser = argparse.ArgumentParser(
        description='Combine PDF pages into a grid.')
    parser.add_argument('input_pdf_path', type=str,
                        help='Path to the input PDF')
    parser.add_argument('output_pdf_path', type=str,
                        help='Path to the output PDF')
    parser.add_argument('--pages_to_merge', type=int, nargs='+', default=None,
                        help='Page numbers to merge, separated by spaces (default: all)')
    parser.add_argument('--grid_shape', type=int, nargs=2, default=[3, 2],
                        help='Number of rows and columns in the grid')
    parser.add_argument('--margins', type=float, nargs=2, default=[0.15, 0.15],
                        help='Horizontal and vertical margins in the grid, as fraction of cell size.')
    args = parser.parse_args()

    merge_pages_into_grid(args.input_pdf_path, args.output_pdf_path,
                          args.pages_to_merge, args.grid_shape, args.margins)
