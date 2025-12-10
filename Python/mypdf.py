from fpdf import FPDF       # pip install fpdf2
from datetime import datetime
# import os
# fpdf.set_global("SYSTEM_TTFONTS", os.path.join(os.path.dirname(__file__),'fonts')
import io
import numbers
from matplotlib.figure import Figure

'''Class PDF derived from fpdf2'''
# https://pyfpdf.readthedocs.io/en/latest/Tutorial/index.html


class PDF(FPDF):
    count_chapter = 0
    count_subchapter = 0
    count_figure = 0
    count_table = 0
    def __init__(self, title, author):
        super().__init__()
        self.title = title
        self.author = author
        self.set_title(title)
        self.set_author(author)
        self.set_margins(24.1, 16.9, 20)    # DIN 5008

        # fontpath = 'C:\\WINDOWS\\FONTS\\arial.ttf'            TODO import of Windows fonts not working
        # self.add_font('Arial', '', fontpath)

    # def header(self):
    #     # Arial bold 15
    #     self.set_font('Arial', 'B', 15)
    #     # Calculate width of title and position
    #     w = self.get_string_width(self.title) + 6
    #     self.set_x((210 - w) / 2)
    #     # Colors of frame, background and text
    #     self.set_draw_color(0, 80, 180)
    #     self.set_fill_color(230, 230, 0)
    #     self.set_text_color(220, 50, 50)
    #     # Thickness of frame (1 mm)
    #     self.set_line_width(1)
    #     # Title
    #     self.cell(w, 9, self.title, 1, 1, 'C', 1)
    #     # Line break
    #     self.ln(10)

    def footer(self):
        # Position at 1.5 cm from bottom
        self.set_y(-15)
        # Helvetica italic 8
        self.set_font('Helvetica', 'I', 8)
        # Text color in gray
        self.set_text_color(128)
        # Document name (left)
        self.cell(0, 10, self.title, 0, 0, 'L')
        # Page number (right)
        self.cell(0, 10, 'Page ' + str(self.page_no()), 0, 0, 'R')

    def introtable(self):
        self.set_font('Helvetica', '', 12)
        creation = datetime.now().strftime("%d.%m.%Y, %H:%M:%S")
        table_data = (('Document', self.title.removesuffix('.pdf')),
                      ('Author', self.author),
                      ('Created', creation))
        # Fill table
        with self.table(first_row_as_headings=False,
                        col_widths=(0.3*self.epw, 0.7*self.epw),
                        align='L') as table:
            for data_row in table_data:
                row = table.row()
                for datum in data_row:
                    row.cell(datum)
        self.ln(10)

    def infotable(self, label, table_data, col_widths=None, precision=2):
        self.set_font('Helvetica', '', 12)
        self.cell(0, 10, label)
        self.ln(10)
        if col_widths is None and len(table_data[0]) == 2:
            col_widths = (0.3 * self.epw, 0.7 * self.epw)
        with self.table(cell_fill_color=None,               # cell_fill_color=None is not working
                        first_row_as_headings=False,
                        col_widths=col_widths,
                        align='L') as table:
            for data_row in table_data:
                row = table.row()
                for datum in data_row:
                    if isinstance(datum, numbers.Number):
                        datum = str(round(datum, precision))        # no problem here
                    row.cell(datum)
        self.ln(10)

    def chapter_title(self, label, shownum=True, num=-1):
        # Helvetica 12
        self.set_font('Helvetica', '', 12)
        # Background color
        self.set_fill_color(200, 220, 255)
        # Title
        if shownum:
            if num < 0:     # default behaviour without specified num
                self.count_chapter += 1
            else:
                self.count_chapter = num
            self.cell(0, 6, 'Chapter %d : %s' % (self.count_chapter, label), 0, 1, 'L', fill=True)
        else:
            self.cell(0, 6, 'Chapter : %s' % label, 0, 1, 'L', fill=True)
        # Reset fill color
        self.set_fill_color(255, 255, 255)
        # Line break
        self.ln(4)

    def chapter_subtitle(self, label, shownum=True, num=-1):
        # Helvetica 12
        self.set_font('Helvetica', 'B', 12)
        # Title
        if shownum:
            if num < 0:
                self.count_subchapter += 1
            else:
                self.count_subchapter = num
            strnum = str(self.count_chapter) + '.' + str(self.count_subchapter)
            self.cell(0, 6, '%s %s' % (strnum, label), 0, 1, 'L')
        else:
            self.cell(0, 6, '%s' % label, 0, 1, 'L')
        # Line break
        self.ln(4)

    # def chapter_body(self, name):
    #     # Read text file
    #     with open(name, 'rb') as fh:
    #         txt = fh.read().decode('latin-1')
    #     # Times 12
    #     self.set_font('Times', '', 12)
    #     # Output justified text
    #     self.multi_cell(0, 5, txt)
    #     # Line break
    #     self.ln()
    #     # Mention in italics
    #     self.set_font('', 'I')
    #     self.cell(0, 5, '(end of excerpt)')

    def chapter_text(self, txt):
        # Helvetica 12
        self.set_font('Helvetica', '', 12)
        # Output justified text
        txt = txt.replace('\n', ' ')
        self.multi_cell(0, 5, txt)
        # Line break
        self.ln()

    def chapter_image(self, fig, label=None, shownum=True, num=-1):
        img = fig.to_image(format='png', scale=2)
        img = io.BytesIO(img)
        self.image(img, x=None, y=None, w=self.epw, h=0, type='', link='')
        if label is not None:
            self.set_font('Helvetica', 'I', 10)
            if shownum:
                if num < 0:
                    self.count_figure += 1
                else:
                    self.count_figure = num
                self.cell(0, 6, 'Figure %d : %s' % (self.count_figure, label), 0, 1, 'C')
            else:
                self.cell(0, 6, 'Figure : %s' % label, 0, 1, 'C')
        self.ln(10)

    def chapter_formula(self, formula, h=20, scale=0.8):
        mm = 1 / 25.4  # millimeters in inches
        fig = Figure(figsize=(self.epw*mm, h*mm))
        gca = fig.gca()
        gca.text(0, 0.5, formula, fontname='Times', fontsize=12, color='black')
        gca.axis("off")

        # Converting Figure to a SVG image:
        img = io.BytesIO()
        fig.savefig(img, format="png")      # TODO svg not working correctly, use .png
        self.image(img, w=scale*self.epw, keep_aspect_ratio=True)

    # def print_chapter(self, num, chaptitle, name):
    #     self.add_page()
    #     self.chapter_title(num, chaptitle)
    #     self.chapter_body(name)
