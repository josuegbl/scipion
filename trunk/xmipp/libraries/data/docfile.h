/***************************************************************************
 *
 * Authors:     Carlos Oscar S. Sorzano (coss@cnb.uam.es)
 *
 * Unidad de  Bioinformatica of Centro Nacional de Biotecnologia , CSIC
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA
 * 02111-1307  USA
 *
 *  All comments concerning this program package may be sent to the
 *  e-mail address 'xmipp@cnb.uam.es'
 ***************************************************************************/

#ifndef DOCFILE_H
#define DOCFILE_H

#include <vector>
#include <string>

#include "matrix1d.h"
#include "matrix2d.h"
#include "selfile.h"

// Forward declaration
class DocFile;

/// @defgroup DocFiles Document Files.

/// @defgroup DocLines Document Lines.
/// @ingroup DocFiles

/** DocFile Line.
 * @ingroup DocLines
 *
 * The Document file is a collection (STL list) of Document lines. This class
 * needn't be accessed in common programs since the DocFile class offers most of
 * the used functions. However, the class is shown in case you may need to
 * access specifically to any of the Document Line functions alone. The DocLines
 * can be either comments or data line, they might be of other internal used
 * types (not assigned or to be discarded), but you may assume that if a line is
 * not a comment nor a data line then it can be skipped.
 */
class DocLine
{
public:
    typedef enum
    {
        NOT_CONSIDERED = -1,
        NOT_ASSIGNED = 0,
        DATALINE = 1,
        COMMENT = 2
    } Line_Type;

private:
    Line_Type line_type;
    string text; ///< if a comment, the whole line
    int key; //< key of this line
    vector< double > data; ///< data of that line

    friend class DocFile; ///< Document file can access anything in this class

public:
    /// @defgroup DocLineConstructors Constructors for Document Lines.
    /// @ingroup DocLines

    /** Empty Constructor.
     * @ingroup DocLineConstructors
     *
     * The document line is created with no type (neither comment or data). You
     * must use the function set_type to assign a type
     */
    DocLine() : line_type(NOT_ASSIGNED)
    {}

    /** Copy constructor.
     * @ingroup DocLineConstructors
     */
    DocLine(const DocLine& line);

    /** Assignment.
     * @ingroup DocLineConstructors
     */
    DocLine& operator=(const DocLine& line);

    /// @defgroup DocLineAccess Component access
    /// @ingroup DocLines

    /** Set existing components.
     * @ingroup DocLineAccess
     *
     * Inside the document line the values are considered as an array (with
     * starting index at 0). With this function you can set any EXISTING
     * position inside the document line. If the position has not been created
     * yet (ie, the vector is smaller) then an exception is thrown.
     *
     * @code
     * line[3] = 3; // line[3] must exist!!
     * @endcode
     */
    double& operator[](int i);

    /** Constant component access.
     * @ingroup DocLineAccess
     *
     * The same as the previous function
     */
    double operator[](int i) const;

    /** Set an existing or not component.
     * @ingroup DocFileAccess
     *
     * If the Document Line is not large enough to hold the required component,
     * then it is resized
     */
    void set(int i, double val);

    /** Set a vector (matrix1D) as Document line.
     * @ingroup DocFileAccess
     *
     * It doesn't matter if it is a row or a column vector. The previous data
     * is overwritten and the key is kept. If it was not a data line, then the
     * new key=0.
     */
    void set(const matrix1D< double >& v);

    /// @defgroup DocLineStructure Structure information.
    /// @ingroup DocLines

    /** Get text of this line.
     * @ingroup DocLineStructure
     *
     * The text is only valid for comment lines
     */
    std::string get_text()
    {
        return text;
    }

    /** Get the key of this line.
     * @ingroup DocLineStructure
     */
    int get_key()
    {
        return key;
    }

    /** Get the number of components.
     * @ingroup DocFileStructure
     *
     * If it is a comment or it hasn't been assigned it returns -1
     */
    int get_no_components()
    {
        if (line_type == DATALINE)
            return data.size();
        else
            return -1;
    }

    /** Empty the document line.
     * @ingroup DocFileStructure
     */
    void clear();

    /** True if current line is a comment.
     * @ingroup DocFileStructure
     */
    int Is_comment()
    {
        return line_type == COMMENT;
    }

    /** True if current line is a comment.
     * @ingroup DocFileStructure
     */
    int Is_data()
    {
        return line_type == DATALINE;
    }

    /** Set type.
     * @ingroup DocFileStructure
     *
     * Only the comment flag is set, the key and possible data are not
     * lost. The comment text is not touched. The valid types are
     * DATALINE, COMMENT, NOT_ASSIGNED, and NOT_CONSIDERED.
     */
    void set_type(Line_Type _line_type)
    {
        line_type = _line_type;
    }

    /** Show a Document Line.
     * TODO document
     */
    friend ostream& operator<<(ostream& o, const DocLine& DL);

    /** Read a Document Line.
     * TODO document
     * An exception is thrown if the line doesn't meet the Document File
     * specifications. First the line is read in a C way, if it fails then the
     * exact Fortran output is tried.
     */
    void read(istream& i);
};

/** Document Files.
 * @ingroup DocFiles
 *
 * The Document Files are the file format for textual information interchange
 * with Spider. These document files are limited to 6 fields if we want a true
 * compatibility but in fact it is not limited and you can use as large document
 * files as you like as long as you keep their use within Xmipp. The files are
 * loaded entirely in memory.
 *
 * The columns in the document file are numbered as 0,1,2,... The first two
 * columns (key and record length) are not taken into account, so column 0 is
 * the first data column.
 *
 * The keys in the document file start at 1 (1,2,3,...) by default but you can
 * change this using the function FirstKey.
 */
class DocFile
{
    FileName fn_doc; ///< Filename
    vector< DocLine > m; ///< All data
    int no_lines; ///< Number of not commented lines
    int first_key; ///< Number of the first file key
    vector< DocLine >::iterator current_line; ///< "pointer" to current line

    // TODO document
    vector< DocLine >::iterator find(int _key);

public:
    /// @defgroup DocFileConstructor DocFile constructors.
    /// @ingroup DocFiles

    /** Empty constructor.
     * @ingroup DocFileConstructor
     */
    DocFile(): fn_doc(""), no_lines(0), first_key(1)
    {
        current_line = m.begin();
    }

    /** Constructor with filename, read from disk.
     * @ingroup DocFileConstructor
     *
     * The given name is loaded as a document file. If it doesn't exist an
     * exception is thrown.
     *
     * @code
     * DocFile DF("angles.doc");
     * @endcode
     */
    DocFile(FileName doc_name)
    {
        first_key = 1;
        read(doc_name);
    }

    /** Copy constructor.
     * @ingroup DocFileConstructor
     *
     * @code
     * DocFile DF2(DF1);
     * @endcode
     */
    DocFile(const DocFile& DF);

    /** Empties the object.
     * @ingroup DocFileConstructor
     *
     * @code
     * DF.clear();
     * @endcode
     */
    void clear();

    /** Reserve for N entries.
     * @ingroup DocFileConstructor
     *
     * It doesn't matter if they are comments or data lines. It's very
     * important to make a reservation if you don't want a memory explosion!!
     *
     * The current line is set to the beginning of the docfile
     */
    void reserve(int N)
    {
        m.reserve(N);
        current_line = m.begin();
    }

    /// @defgroup DocFileOperator DocFile operators.
    /// @ingroup DocFiles

    /** Assignment
     * @ingroup DocFileOperator
     *
     * @code
     * DF2 = DF1;
     * @endcode
     */
    DocFile& operator=(const DocFile &DF);

    /** Another function for assigment.
     * @ingroup DocFileOperator
     */
    void assign (const DocFile &DF);

    /** Assignment from matrix.
     * @ingroup DocFileOperator
     *
     * The old information in the document file is lost.
     *
     * @code
     * matrix2D< double > A(30, 3);
     * A.init_random(0, 180);
     * DF = A;
     * @endcode
     */
    DocFile& operator=(const matrix2D< double >& A);

    /** Another function for assignment from matrix.
     * @ingroup DocFileCOprator
     */
    void assign(const matrix2D< double >& A);

    /** Show a document file.
     * @ingroup DocFileOperator
     *
     * A new line is printed at the end.
     *
     * @code
     * std::cout << DF;
     * @endcode
     */
    friend ostream& operator<<(ostream& o, const DocFile& DF);

    /** Show a given line.
     * @ingroup DocFileOperator
     *
     * This function shows you the line with the given key, if no key is given
     * then the current line is shown. After showing the line a newline
     * character is shown.
     *
     * @code
     * DF.show_line(); // Show current line
     * DF.show_line(3);  // Show line with key=3
     * @endcode
     */
    void show_line(ostream& o, int key=-1);

    /** Show everything in a document file.
     * @ingroup DocFileOperator
     */
    void debug();

    /// @defgroup DocFileDisk Managing files in disk.
    /// @ingroup DocFiles

    /** Read a file from disk.
     * @ingroup DocFileDisk
     *
     * The old information on the variable is overwritten. An exception is
     * thrown if the file doesn't exist. Lines which do not fit the file format
     * are ignored.
     *
     * First the file is tried to be read in a loose C way. If it fails then the
     * Fortran style is tried. It is especially useful for reading files coming
     * from Spider because there might be no space between two fields.
     *
     * @code
     * DF.read("angles.doc");
     * @endcode
     */
    void read(FileName _name, int overrinding=1);

    /** Append a file from disk to an already read one.
     * @ingroup DocFileDisk
     *
     * The old information on the variable is not lost. All lines in the
     * document file to be read are appened at the end of the already read one.
     *
     * @code
     * DF.read("angles1.doc");
     * DF.append("angles2.doc");
     * @endcode
     */
    void append(FileName _name)
    {
        read(_name,0);
    }

    /** Write a document file to disk.
     * @ingroup DocFileDisk
     *
     * If you give a name then it becomes like a "Save as ..." and from this
     * point on the name of the document file has changed. The file keys are
     * renumerated before being saved.
     *
     * @code
     * DF.write(); // Save
     * DF.write("angles3.doc"); // Save as ...
     * @endcode
     */
    void write(FileName _name="");

    /// @defgroup DocFilePointer Moving the current line pointer.
    /// @ingroup DocFiles

    /** Go to the beginning of the file.
     * @ingroup DocFilePointer
     *
     * Moves the pointer to the first line of the file either it is a comment or
     * a data line.
     *
     * @code
     * DF.go_beginning();
     * @endcode
     */
    void go_beginning()
    {
        current_line = m.begin();
    }

    /** Go to the first data line.
     * @ingroup DocFilePointer
     *
     * Moves the pointer to the first data line in the file.
     *
     * @code
     * DF.go_first_data_line();
     * @endcode
     */
    void go_first_data_line()
    {
        go_beginning();
        adjust_to_data_line();
    }

    /** Adjust pointer to a data line.
     * @ingroup DocFilePointer
     *
     * If the current line is a data line, nothing is done, else the current
     * line pointer is moved until it is pointing a data line.
     *
     * @code
     * DF.adjust_to_data_line();
     * @endcode
     */
    void adjust_to_data_line();

    /** Moves pointer to next line.
     * @ingroup DocFilePointer
     *
     * It doesn't matter if next line is a comment or a data line.
     *
     * @code
     * DF.next();
     * @endcode
     */
    void next()
    {
        if (current_line != m.end())
            current_line++;
    }

    /** Moves pointer to previous line.
     * @ingroup DocFilePointer
     *
     * It doesn't matter if previous line is a comment or a data line.
     *
     * @code
     * DF.previous();
     * @endcode
     */
    void previous()
    {
        if (current_line != m.begin())
            current_line--;
    }

    /** Moves pointer to next data line.
     * @ingroup DocFilePointer
     *
     *
     * @code
     * DF.next_data_line();
     * @endcode
     */
    void next_data_line()
    {
        jump(1);
    }

    /** Jump over a number of data lines.
     * @ingroup DocFilePointer
     *
     * Starting from the current_line "pointer" this function skips a given
     * number of entries. For instance, jump over 1 data line is to jump to the
     * next data line. Jump over 2 data lines is to jump to the next of the next
     * data line, and so on. The number of images to jump must always be
     * positive, the jumpings cannot be done backwards.
     *
     * @code
     * DF.jump(2); // Jump over 2 data lines
     * @endcode
     */
    void jump(int how_many);

    /** Move "pointer" to a certain line.
     * @ingroup DocFilePointer
     *
     * This function searches for the line with the given key, and locates the
     * current line "pointer" pointing to that line. If the key is not present
     * in the document file, then the pointer is moved to the end of the
     * document file. You can check this situation using eof().
     *
     * @code
     * DF.search(700);
     * @endcode
     */
    void search(int _key)
    {
        current_line = find(_key);
    }

    /** Search the entire file for the given comment.
     * @ingroup DocFilePointer
     *
     * If found, place the pointer to the next data line and return 1.
     * Otherwise, return 0.
     *
     * @code
     * if (DF.search_comment(fn_img))
     *     rot=DF(0);
     * @endcode
     */
    int search_comment(std::string comment);

    /** Extract the comment's images' SelFile.
     * @ingroup DocFilePointer
     *
     * For NewXmipp-type Docfiles only
     *
     * @code
     * DF.get_selfile(SF);
     * @endcode
     */
    void get_selfile(SelFile& SF);

    /** Move "pointer" to a certain line.
     * @ingroup DocFilePointer
     *
     * This function searches for the line with the given key, and locate the
     * current line "pointer" pointing to that line. If the key is not present
     * in the document file, then the pointer is placed pointing to the line
     * with the smaller key greater than the given one. For instance, if you ask
     * for key 700, but there is no key from 600-800 (inclusive), then the
     * current_line is pointing to line number 801.
     *
     * @code
     * DF.locate(700);
     * @endcode
     */
    void locate(int _key);

    /** True if current line "pointer" is at the end of file.
     * @ingroup DocFilePointer
     *
     * @code
     * if (SF.eof())
     *     cout << "The document file is over\n";
     * @endcode
     */
    int eof()
    {
        return current_line == m.end();
    }

    /// @defgroup DocFileInfo Getting information.
    /// @ingroup DocFiles

    /** Returns the name of the file
     * @ingroup DocFileInfo
     */
    std::string name() const
    {
        return fn_doc;
    }

    /** True if the key is inside the document file.
     * @ingroup DocFileInfo
     *
     * The current line "pointer" is not modified.
     *
     * @code
     * if (DF.exists(700))
     *     std::cout << "key 700 exists in the document file" << std::endl;
     * @endcode
     */
    int exists(int _key)
    {
        return find(_key) != m.end();
    }

    /** Number of columns of the first data line.
     * @ingroup DocFileInfo
     *
     * The current line "pointer" is not modified.
     */
    int FirstLine_ColNo();

    /** Number of data lines inside a document file.
     * @ingroup DocFileInfo
     *
     * @code
     * std::cout << "There are " << DF.dataLineNo();
     * std::cout << " data lines" << std::endl;
     * @endcode
     */
    int dataLineNo()
    {
        return no_lines;
    }

    /** Returns the number of lines within a file.
     * @ingroup DocFileInfo
     *
     * This function gives the total number of lines (including comments)
     * within a file.
     *
     * @code
     * std::cout << "There are " << DF.LineNo();
     * std::cout << " lines in this file" << std::endl;
     * @endcode
     */
    int LineNo()
    {
        return m.size();
    }

    /** Get last key of the current file.
     * @ingroup DocFileInfo
     *
     * If the file is empty then returns 0.
     *
     * @code
     * last_key = DF.get_last_key();
     * @endcode
     */
    int get_last_key();

    /** Get the key of the current line.
     * @ingroup DocFileInfo
     *
     * If the current line "pointer" is at the end of the file or is pointing to
     * a comment then 0 is returned.
     *
     * @code
     * current_key = DF.get_current_key();
     * @endcode
     */
    int get_current_key()
    {
        return (*current_line).key;
    }

    /** Get first key of the file.
     * @ingroup DocFileInfo
     *
     * This is not the first existing key of the file (the line with this key
     * might be deleted), but the key where all renumerations start.
     *
     * @code
     * std::cout << "First key : " << DF.FirstKey() << std::endl;
     * @endcode
     */
    int FirstKey() const
    {
        return first_key;
    }

    /** Set first key of the file.
     * @ingroup DocFileInfo
     *
     * This is not the first existing key of the file (the line with this key
     * might be deleted), but the key where all renumerations start. No
     * renumeration is performed at this point.
     *
     * @code
     * DF.FirstKey()=700;
     * @endcode
     */
    int& FirstKey()
    {
        return first_key;
    }

    /** Another function for set first key of the file.
     * @ingroup DocFileInfo
     */
    void set_FirstKey(int& _first_key)
    {
        first_key = _first_key;
    }

    /** Get the number of the values in the current line.
     * @ingroup DocFileInfo
     *
     * If the current line "pointer" is at the end of the file or is pointing to
     * a comment then 0 is returned.
     *
     * @code
     * valNo = DF.get_current_valNo();
     * @endcode
     */
    int get_current_valNo()
    {
        return (*current_line).data.size();
    }

    /** Constant access to a value in the current line.
     * @ingroup DocFileInfo
     *
     * This function allows you access to the values inside the current line of
     * the document file. The column numbers start at 0.
     *
     * @code
     * std::cout << "Value at first column of the current line = ";
     * std::cout << DF(0) << std::endl;
     * @endcode
     */
    double operator()(int i) const
    {
        return (*current_line)[i];
    }

    /** Constant access to a value in the current line.
     * @ingroup DocFileInfo
     *
     * This function allows you access to the values inside any line of the
     * document file. The column numbers start at 0. The current line "pointer"
     * is not modified. If the data line hasn't got already enough space to hold
     * the desired column, space is allocated automatically for that line. Be
     * careful that the not assigned columns contain garbage and should be
     * assigned sooner or later.
     *
     * @code
     * DF(700, 0) = 4;
     * @endcode
     */
    double operator()(int _key, int i);

    /** Get angles on key i.
     * @ingroup DocFileInfo
     *
     * You must specify the order in which they are written in the file with the
     * labels "rot", "tilt" and "psi"
     */
    void get_angles(int _key, double& rot, double& tilt, double& psi,
                    const string& ang1, const string& ang2, const string& ang3);

    /** Get angles on key i (second triad of Euler angles).
     * @ingroup DocFileInfo
     *
     * You must specify the order in which they are written in the file with the
     * labels "rot", "tilt" and "psi"
     */
    void get_angles1(int _key, double& rot, double& tilt, double& psi,
                     const string& ang1, const string& ang2,
                     const string& ang3);

    /** Get angles on key i (third triad of Euler angles)
     * @ingroup DocFileInfo
     *
     * You must specify the order in which they are written in the file with the
     * labels "rot", "tilt" and "psi"
     */
    void get_angles2(int _key, double& rot, double& tilt, double& psi,
                     const string& ang1, const string& ang2,
                     const string& ang3);

    /** Set angles on key i.
     * @ingroup DocFileInfo
     *
     * You must specify the order in which they are written in the file with the
     * labels "rot", "tilt" and "psi"
     */
    void set_angles(int _key, double rot, double tilt, double psi, const string&
                    ang1, const string& ang2, const string& ang3);

    /** Set a value in the current line.
     * @ingroup DocFileInfo
     *
     * This function allows you to set values inside the current line of the
     * document file. The column numbers start at 0. If the data line hasn't got
     * already enough space to hold the desired column, space is
     * allocated automatically for that line. Be careful that the not assigned
     * columns contain 0's and should be assigned sooner or later. The following
     * example sets a 4 at column 0 of the current line. If the current line is
     * pointing to the end of the document file then a new line without key is
     * added with the asked value at the desired column. The current line is
     * moved again to the end of the file.
     *
     * @code
     * DF.set(0, 4);
     * @endcode
     */
    void set(int i, double val);

    /** Set a value in the given line.
     * @ingroup DocFileInfo
     *
     * This function allows you access to the values inside any line of the
     * document file (the line with the key must exist). The column numbers
     * start at 0. The current line "pointer" is not modified. If the data line
     * hasn't got already enough space to hold the desired column, space is
     * allocated automatically for that line. Be careful that the not assigned
     * columns contain 0's and should be assigned sooner or later. The following
     * example sets a 4 at column 0 of the line with key 700.
     *
     * @code
     * DF.set(700, 0, 4);
     * @endcode
     */
    void set(int _key, int i, double val);

    /** Returns current line as a Document line.
     * @ingroup DocFileInfo
     */
    DocLine get_current_line()
    {
        return *current_line;
    }

    /// @defgroup DocFileModify Modifying the document file.
    /// @ingroup DocFiles

    /** Renumerate keys
     * @ingroup DocFileModify
     *
     * The keys are renumerated starting from the file first_key (by default, 1)
     * at the beginning of the file. The current line "pointer" is not modified
     *
     * @code
     * DF.renum();
     * @endcode
     */
    void renum();

    /** Removes a line or several lines from the document file.
     * @ingroup DocFileModify
     *
     * This function searches for a key (or a range of keys) in the document
     * file, if it is found then the corresponding line(s) is(are) deleted. If
     * the line is actually being pointed by the current line, then the current
     * line is now the following one. Keys are NOT renumerated at the end, in
     * this way the keys remain the same as before and you know, for instance in
     * the following examples, that keys 706- are still valid.
     *
     * @code
     * DF.remove(700);
     * DF.remove(700,705); // both limits are removed too
     * @endcode
     */
    void remove(int _key0, int _keyF=-1);

    /** Removes actual line.
     * @ingroup DocFileModify
     *
     * This function removes the current line, either it is a comment or a data
     * line. The current line "pointer" is moved to the following line in the
     * file. Keys are NOT renumerated at the end (see previous function for an
     * explanation for this design option.
     *
     * @code
     * DF.remove_current();
     * @endcode
     */
    void remove_current();

    /** Insert empty data lines before current line.
     * @ingroup DocFileModify
     *
     * The current line is still pointing to the same line as it was before
     * entering the function. The first key assigned to the new lines is
     * returned. After the insertion the file keys are renumerated. By default
     * only one line is inserted but you can tell to insert more lines.
     *
     * @code
     * new_key = DF.insert_data_line();
     * @endcode
     */
    int insert_data_line(int no_lines_to_insert=1);

    /** Insert a data line before current line.
     * @ingroup DocFileModify
     *
     * The current line is still pointing to the same line as it was before
     * entering the function. The key assigned to the new line is returned.
     * After the insertion the file keys are renumerated. It doesn't matter if
     * the vector initially is a row or a column, it will always be inserted as
     * a row vector, ie, as a line in the document file. Notice that several
     * succesive insertions are "queued" one after another, even if the document
     * file is empty at the beginning.
     *
     * @code
     * // after many insertions the document file stage would be
     *  line             n-1
     *  insertion         1
     *  insertion         2
     *  ...
     *  insertion         m
     *  line              n
     *
     * matrix1D< double > proj_angles(3);
     * ...
     * new_key = DF.insert_data_line(proj_angles);
     * @encode
     */
    int insert_data_line(const matrix1D< double >& v);

    /** Insert a comment before the current line.
     * @ingroup DocFileModify
     *
     * Comments must not start with any special character since a ";" is
     * automatically added at the beginning of the line. The current line is
     * still pointing to the same line as it was before entering the function.
     * After the insertion the file keys are NOT renumerated.
     *
     * @code
     * DF.insert_comment("This is a comment");
     * @endcode
     */
    void insert_comment(string comment);

    /** Insert a document line before the current line.
     * @ingroup DocFileModify
     *
     * The document line is inserted in the file as it is, the key is modified
     * by a renumeration (if the document line is a comment then no renumeration
     * is performed).
     *
     * The current line is still pointing to the same line as it was before
     * entering the function. After the insertion the file keys are renumerated
     * if the Document Line is a data line (then the assigned key is returned),
     * if the inserted line is a comment then 0 is returned.
     *
     * @code
     * DF.insert_comment(DF2.get_current_line());
     * @endcode
     */
    int insert_line(const DocLine& DL);

    /** Append empty data lines before at the end of the file.
     * @ingroup DocFileModify
     *
     * The current line is still pointing to the same line as it was before
     * entering the function. The first key assigned to the new lines is
     * returned. After the appending no renumeration is needed. By default only
     * one line is appended but you can tell to append more lines.
     *
     * @code
     * new_key = DF.append_data_line();
     * @endcode
     */
    int append_data_line(int no_lines_to_append=1);

    /** Append a data line at the end of the file.
     * @ingroup DocFileModify
     *
     * This function is exactly the same as the insertion one, but no
     * renumeration is performed at the end, so some speed is gained. The key
     * assigned goes on being returned and the current_line pointer is not
     * moved, neither.
     *
     * @code
     * matrix1D< double > proj_angles(3);
     * ...
     * new_key = DF.append_data_line(proj_angles);
     * @endcode
     */
    int append_data_line(const matrix1D< double >& v);

    /** Append angles.
     * @ingroup DocFileModify
     *
     * You must specify the order with "rot","tilt", or "psi"
     */
    int append_angles(double rot, double tilt, double psi, const string& ang1,
                      const string& ang2, const string& ang3);

    /** Append angles, using 2 triads of Euler angles.
     * @ingroup DocFileModify
     *
     * You must specify the order with "rot","tilt", or "psi"
     */
    int append_angles(double rot,  double tilt,  double psi, double rot1, double
                      tilt1, double psi1, const string& ang1, const string&
                      ang2, const string& ang3);

    /** Append angles, using three triads of Euler angles.
     * @ingroup DocFileModify
     *
     * You must specify the order with "rot","tilt", or "psi"
     */
    int append_angles(double rot, double tilt, double psi, double rot1,
                      double tilt1, double psi1, double rot2, double tilt2,
                      double psi2, const string& ang1, const string& ang2,
                      const string& ang3);

    /** Append a comment at the end of the file.
     * @ingroup DocFileModify
     *
     * This function is exactly the same as the insertion one, but this time the
     * line is added at the end. The current line pointer is not moved and no
     * renumeration is performed.
     *
     * @code
     * DF.append_comment("This is a comment");
     * @endcode
     */
    void append_comment(const string& comment);

    /** Append a document line at the end of the file.
     * @ingroup DocFileModify
     *
     * This function is exactly the same as the insertion one, but no
     * renumeration is performed at the end, so some speed is gained. The key
     * assigned goes on being returned (if the appended line is a data line, if
     * not 0 is returned) and the current_line pointer is not moved, neither.
     *
     * @code
     * new_key = DF.append_line(DF2.get_current_line());
     * @endcode
     */
    int append_line(DocLine& DL);

    /** Deletes all comments from the document file.
     * @ingroup DocFileModify
     *
     * The current line "pointer" is moved to the beginning of the file.
     *
     * @code
     * DF.clean_comments();
     * @endcode
     */
    void clean_comments();

    /// @defgroup DocFileHelpul Helpful functions.
    /// @ingroup DocFiles

    /** Alter order in document file.
     * @ingroup DocFileHelpful
     *
     * A new document file with all the lines of the actual object is created,
     * but this time all data lines are resorted in a random order. The comments
     * of the original document file are lost in the new copy. The current line
     * of the resulting document file is placed at the beginning of the file and
     * the keys are renumerated
     *
     * @code
     * DF2 = DF1.randomize();
     * @endcode
     */
    DocFile randomize();

    /** Discard randomly N lines.
     * @ingroup DocFileHelpful
     *
     * A set of N data lines are REMOVED from the actual document file. If N is
     * equal or greater than the actual number of lines within the file, all
     * lines are discarded. Comments are kept at their original positions. The
     * current line of the resulting document file is placed at the beginning of
     * the file.
     *
     * @code
     * DF2 = DF1.random_discard(3);
     * @endcode
     */
    DocFile random_discard(int N);

    /** Column to vector.
     * @ingroup DocFileHelpful
     *
     * This function produces a double matrix1D which is composed by all the
     * components of a certain column inside the document file. If a given line
     * hasn't got enough data to fill that column (for instance, the line is 3
     * values long and we are asking for column 3 (remember that column
     * numbering start at 0)) then those values are supplied by the function as
     * 0's.
     *
     * @code
     * matrix1D< double > tilt_angle = DF.column(1);
     * @endcode
     */
    matrix1D< double > col(int _col);

    /** Row to vector.
     * @ingroup DocFileHelpful
     *
     * This function produces a double matrix1D which is composed by all the
     * components of a certain line inside the document file. If the key doesn't
     * exist then an empty vector is returned.
     *
     * @code
     * matrix1D< double > tilt_angle = DF.row(700);
     * @endcode
     */
    matrix1D< double > row(int _key);

    /** Vector to column.
     * @ingroup DocFileHelpful
     *
     * This function sets all values in the given column as the values given in
     * the vector. If a data line hasn't got enough space to hold that column,
     * 0's are added. If there are less values in the vector than data lines in
     * the document file an exception is thrown but the values are set. If there
     * are more values in the vector than data lines in the document file, new
     * data lines are added at the end. A renumeration is performed at the end
     * of the setting, and the current line "pointer" is moved to the beginning
     * of the file.
     *
     * @code
     * DF.set_column(1, tilt_angle);
     * @endcode
     */
    void setCol(int _col, matrix1D< double >& v);

    /** Apply a function to all data lines within a key range.
     * @ingroup DocFileHelpful
     *
     * The function must take as input a row double matrix1D and produce another
     * row double matrix1D which will be assigned instead of the input data
     * line. There is a renumeration of lines at the end and the current line
     * "pointer" is moved to the beginning of the file. The document file must
     * entered in the routine previously renumerated.
     *
     * @code
     * // This function enlarges the row vector adding 360 at the beginning
     * void enlarge(const matrix1D< double >&in, matrix1D< double >& out)
     * {
     *     out.resize(in.get_dim() + 1);
     *     out(0) = 360;
     *     for (int i=1; i<out.get_dim(); i++)
     *         out(i)=in(i-1);
     * }
     *
     * ...
     * DF.for_all_lines(&enlarge); // Apply to all lines
     *
     * // Apply to lines with key greater or equal to 30
     * DF.for_all_lines(&enlarge, 30);
     *
     * // Apply to all lines with 30 <= key <= 40
     * DF.for_all_lines(&enlarge, 30, 40);
     *
     * @endcode
     */
    void for_all_lines(void (*f) (const matrix1D< double >&, matrix1D< double >&
                                 ), int key0=-1, int keyF=-1);

    /** Apply a function to a data line with a given key.
     * @ingroup DocFileHelpful
     *
     * The function requirements are the same as for the preceeding function
     * (see there for a more detailed explanation). If the key doesn't exist
     * nothing is done. It is based on for_all_lines, so there is a renumeration
     * at the end.
     *
     * @code
     * DF.for_one_line(&function1, 30);
     * // Apply to line with key=30
     * @endcode
     */
    void for_one_line(void (*f) (const matrix1D< double >&, matrix1D< double >&)
                      , int key0)
    {
        for_all_lines(f, key0, key0);
    }


    /** Apply a function to a column between a given key range.
     * @ingroup DocFileHelpful
     *
     * The function must take a double and return a double. The function is
     * applied to the given column (if nothing is said then it is applied to all
     * columns) between the given range of keys. There is a renumeration at the
     * end, and the current line "pointer" is moved to the beginning of the
     * file.
     *
     * The document file must entered in the routine previously renumerated.
     *
     * @code
     * double adjust_angle(double ang) { return PI/2-ang; }
     *
     * ...
     * DF.for_column(&adjust_angle); // Apply to all values in all columns
     * DF.for_column(&adjust_angle, 1); // Apply to all values in column 1
     *
     * // Apply to all values in column 1 with key greater or equal to 30
     * DF.for_column(&adjust_angle, 1, 30);
     *
     * // Apply to all values in column 1 for which 30 <= key <= 40
     * DF.for_column(&adjust_angle, 1, 30, 40);
     * @endcode
     */
    void for_column(double (*f) (double), int _col=-1, int _key0=-1, int
                    _keyF=-1);
};

/** Read Document File with Euler angles.
 * @ingroup DocFileHelpful
 *
 * This function reads a document file with Euler angles. The Euler angles can
 * be in any order in the file (the order must be specified in the function
 * call, see below) but the result is always reordered to (rot, tilt, psi). The
 * function returns the number of total data lines inside the document file.
 *
 * The order specification is performed using the labels "rot", "tilt" and
 * "psi". You must specify a label for each column. The following example reads
 * "angle.doc" into DF. The order of columns inside angle.doc is (psi, tilt,
 * rot) but DF will be (rot,tilt,psi).
 *
 * @code
 * std::cout << "There are ";
 * std::cout << read_Euler_document_file("angle.doc","psi","tilt","rot",DF);
 * std::cout << " sets of angles" << std::endl
 * @endcode
 */
int read_Euler_document_file(FileName fn, string ang1, string ang2, string ang3,
                             DocFile& DF);

/** Select images from a selfile meeting some condition.
 * @ingroup DocFileHelpful
 *
 * If an image is discarded in the selfile it remains discarded in spite of it
 * meets the condition. If it is active then it is discarded if it doesn't meet
 * the condition. Conditions are set on the column col of the DocFile.
 * Conditions are specified by en_limit (enable limit) and the corresponding
 * value.
 *
 * It is not checked that the number of images (active and discarded) and the
 * number of data linea in the DocFile match.
 */
void select_images(DocFile& DF, SelFile& SF, int col, bool en_limit0, double
                   limit0, bool en_limitF, double limitF);

#endif
