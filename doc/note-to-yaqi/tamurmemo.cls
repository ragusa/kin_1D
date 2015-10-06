%%---------------------------------------------------------------------------%%
%% RMEMO.CLS:
%% ----------
%% 
%%---------------------------------------------------------------------------%%
\NeedsTeXFormat{LaTeX2e}
\ProvidesClass{tamurmemo}
        [07/12/05 TAMU Research Memo Class, J. Morel]
\DeclareOption*{\PassOptionsToClass{\CurrentOption}{letter}}
\ProcessOptions
\LoadClass{letter}
\RequirePackage{graphicx}
\usepackage[dvips]{color}
\definecolor{maroon}{rgb}{0.5,0,0}
%% Modify the normal text area on the page.
\setlength{\textwidth}{6.5in}
\setlength{\textheight}{8.0in}
\setlength{\oddsidemargin}{0in}
\setlength{\topmargin}{0.65in}
\setlength{\headheight}{3.5ex}
\setlength{\headsep}{18pt}
%% Add new commands for the LANL memo class.
\newcommand*{\collegename}[1]{\def\colname{#1}}
\newcommand*{\deptname}[1]{\def\deptname{#1}}
\newcommand*{\fromms}[1]{\def\fromwho{#1}}
\newcommand*{\originator}[1]{\def\originatr{\uppercase{#1}}}
\newcommand*{\typist}[1]{\def\typst{#1}}
\newcommand*{\toms}[1]{\def\towho{#1}\gettoname{#1}}
\newcommand*{\thru}[1]{\def\through{#1}}
\newcommand*{\refno}[1]{\def\refnum{#1}}
\newcommand*{\subject}[1]{\def\subjectline{#1}}
\newcommand{\distribution}[1]{\def\distrib{#1}}
\newcommand{\cy}[1]{\def\cyy{#1}}
\newcommand{\enc}[1]{\def\encl{#1}}
\newcommand{\encas}{\def\encass{Enc. a/s}}
\newcommand{\attachments}[1]{\def\attach{#1}}
\newcommand{\attachmentas}{\def\attachas{Attachment a/s}}
\newcommand{\attachmentsas}{\def\attachsas{Attachments a/s}}
% Revise opening command to add logo and info block
\renewcommand*{\opening}{
\thispagestyle{firstpage}
$ $\\
\begin{center}
\vspace{-2.0in}
%%{\includegraphics[width=0.75in]{psfiles/TAMUlogo.eps}}\\
%%%{\includegraphics[width=0.75in]{psfiles/uniseal_1color.eps}}\\
%{\includegraphics[width=0.1in]{seal.png}}\\
%\vspace{0.2in}
%{\rmfamily\bfseries\fontsize{16}{16}\selectfont \rule{0pt}{21pt}TEXAS A\&M UNIVERSITY}\\
\textcolor{maroon}{
{\large \bf TEXAS A\&M UNIVERSITY}\\
{\large  \colname}\\
{\large \deptname}\\
}
{\sffamily\bfseries\fontsize{18}{18}\selectfont \rule{0pt}{19pt}Research \kern.02em Memo}
\end{center}
\par
\begin{minipage}[t]{4.0in}
\begin{tabular}[t]{ll}
{\scriptsize From:}  & \fromwho\\
  \ifx\through\@empty
  \else
  {\scriptsize Thru:} & \through\\
  \fi
  \ifx\refnum\@empty
  \else
  {\scriptsize Symbol:}&  \refnum
  \fi
\end{tabular}
\end{minipage}
\begin{minipage}[t]{2.5in}
\hfill 
\begin{tabular}[t]{ll}
{\scriptsize To:} & \towho\\  
{\scriptsize Date:} & \@date
\end{tabular}
\end{minipage}
\par
\vspace{24pt}
\textbf{\large Subject: \subjectline}
\par
}
% Revise closing command to add additional memo options.
\renewcommand*{\closing}{
   \stopbreaks
   \par
   \noindent
   \originatr:\typst
   \ifx\attach\@empty\else
     \par 
     \noindent
     Attachments:\\
     \attach
   \fi
   \ifx\attachas\@empty\else
     \par 
     \noindent
     \attachas
   \fi
   \ifx\attachsas\@empty\else
     \par 
     \noindent
     \attachsas
   \fi
   \ifx\encl\@empty\else
     \par 
     \noindent
     Enclosures:\\
     \encl
   \fi
   \ifx\encass\@empty\else
     \par 
     \noindent
     \encass
   \fi
   \ifx\distrib\@empty\else
     \par 
     \noindent
     Distribution:\\
     \distrib
   \fi
   \ifx\cyy\@empty\else
     \par 
     \noindent
     Cy:\\
     \cyy
   \fi
   }
% Procedure to strip the toname out of "toname, group, MS".
\def\Dist{Distribution}
\def\gettoname#1{\ifx\towho\Dist\let\toname=\towho\else\strip #1\par\fi}
%\long\def\strip #1,#2@@@{\def\toname{#1}\def\tojunk{#2}}
\long\def\strip #1,#2\par{\def\toname{#1}\def\tojunk{#2}}
% Revise Headings
% \refnum
\if@twoside
  \def\ps@headings{%
      \let\@oddfoot\@empty\let\@evenfoot\@empty
      \def\@oddhead{\slshape%
        \begin{tabular}{l}
          \headtoname{} \ignorespaces\toname \\
          {footnotesize\refnum}
        \end{tabular}
                    \hfil --\thepage--
                    \hfil \@date}%
      \let\@evenhead\@oddhead}
\else
  \def\ps@headings{%
      \let\@oddfoot\@empty
      \def\@oddhead{\slshape%
        \begin{tabular}{l}
          \headtoname{} \ignorespaces\toname \\
          {\footnotesize\refnum}
        \end{tabular}
                    \hfil --\thepage--
                    \hfil \@date}}
\fi
% Initializations
\distribution{}
\cy{}
\enc{}
\attachments{}
\def\attachas{}
\def\attachsas{}
\def\encass{}
\def\toname{}
\def\tojunk{}
\def\through{}
\def\refnum{}
\pagestyle{headings}
%--Bibliographystyle
\newdimen\bibindent
\setlength\bibindent{1.5em}
\newcommand\refname{References}
\newenvironment{thebibliography}[1]
%     {\section*{\refname
%        \@mkboth{\MakeUppercase\refname}{\MakeUppercase\refname}}%
      {\normalfont\bfseries\refname\normalfont
      \list{\@biblabel{\@arabic\c@enumiv}}%
           {\settowidth\labelwidth{\@biblabel{#1}}%
            \leftmargin\labelwidth
            \advance\leftmargin\labelsep
            \@openbib@code
            \usecounter{enumiv}%
            \let\p@enumiv\@empty
            \renewcommand\theenumiv{\@arabic\c@enumiv}}%
      \sloppy\clubpenalty4000\widowpenalty4000%
      \sfcode`\.\@m}
     {\def\@noitemerr
       {\@latex@warning{Empty `thebibliography' environment}}%
      \endlist}
\newcommand\newblock{\hskip .11em\@plus.33em\@minus.07em}
\let\@openbib@code\@empty
%--Figures and Tables
\newcommand\figurename{FIG.}
\newcommand\tablename{TABLE}
\newcounter{figure}
\renewcommand\thefigure{\@arabic\c@figure}
\def\fps@figure{tbp}
\def\ftype@figure{1}
\def\ext@figure{lof}
\def\fnum@figure{\figurename~\thefigure}
\newenvironment{figure}
               {\@float{figure}}
               {\end@float}
\newenvironment{figure*}
               {\@dblfloat{figure}}
               {\end@dblfloat}
\newcounter{table}
\renewcommand\thetable{\@arabic\c@table}
\def\fps@table{tbp}
\def\ftype@table{2}
\def\ext@table{lot}
\def\fnum@table{\tablename~\thetable}
\newenvironment{table}
               {\@float{table}}
               {\end@float}
\newenvironment{table*}
               {\@dblfloat{table}}
               {\end@dblfloat}
\newlength\abovecaptionskip
\newlength\belowcaptionskip
\setlength\abovecaptionskip{10\p@}
\setlength\belowcaptionskip{0\p@}
\long\def\@makecaption#1#2{%
  \vskip\abovecaptionskip
  \sbox\@tempboxa{#1: #2}%
  \ifdim \wd\@tempboxa >\hsize
    #1: #2\par
  \else
    \global \@minipagefalse
    \hb@xt@\hsize{\hfil\box\@tempboxa\hfil}%
  \fi
  \vskip\belowcaptionskip}
%--Sections
%-DEFINITION OF SECTIONS
\def\@startsection#1#2#3#4#5#6{%
 \if@noskipsec \leavevmode \fi
 \par \@tempskipa #4\relax
 \@afterindenttrue
 \ifdim \@tempskipa <\z@ \@tempskipa -\@tempskipa \@afterindentfalse\fi
 \if@nobreak \everypar{}\else
     \addpenalty\@secpenalty\addvspace\@tempskipa\fi
 \@ifstar{\@dblarg{\@sect{#1}{\@m}{#3}{#4}{#5}{#6}}}%
         {\@dblarg{\@sect{#1}{#2}{#3}{#4}{#5}{#6}}}%
}
\def\@sect#1#2#3#4#5#6[#7]#8{%
  \edef\@toclevel{\ifnum#2=\@m 0\else\number#2\fi}%
  \ifnum #2>\c@secnumdepth \let\@secnumber\@empty
  \else \@xp\let\@xp\@secnumber\csname the#1\endcsname\fi
 \ifnum #2>\c@secnumdepth
   \let\@svsec\@empty
 \else
    \refstepcounter{#1}%
    \edef\@svsec{\ifnum#2<\@m
       \@ifundefined{#1name}{}{%
         \ignorespaces\csname #1name\endcsname\space}\fi
       \@nx\textup{%
      \@nx\mdseries
         \csname the#1\endcsname.}\enspace
    }%
  \fi
  \@tempskipa #5\relax
  \ifdim \@tempskipa>\z@ % then this is not a run-in section heading
    \begingroup #6\relax
    \@hangfrom{\hskip #3\relax\@svsec}{\interlinepenalty\@M #8\par}%
    \endgroup
    %\ifnum#2>\@m \else \@tocwrite{#1}{#8}\fi
  \else
  \def\@svsechd{#6\hskip #3\@svsec
    \@ifnotempty{#8}{\ignorespaces#8\unskip
       \@addpunct.}%
    %\ifnum#2>\@m \else \@tocwrite{#1}{#8}\fi
  }%
  \fi
  \global\@nobreaktrue
  \@xsect{#5}}
\let\@ssect\relax
%-DEFINITION OF SECTIONS
\setcounter{secnumdepth}{3} \newcounter{section}
\newcommand{\section}{\@startsection {section} {1} {0em}
  {\baselineskip} {\baselineskip} {\center\normalfont\large\bfseries}}
\newcounter{subsection}[section]
\newcommand{\subsection}{\@startsection {subsection} {2} {0em}
  {.5\baselineskip} {.5\baselineskip}
  {\normalfont\bfseries\normalsize}}
\newcounter{subsubsection}[subsection]
\newcommand{\subsubsection}{\@startsection {subsubsection} {3} {0em}
  {-.5\baselineskip} {-2ex} {\normalfont\itshape\normalsize}}
\renewcommand\thesection {\@arabic\c@section}
\renewcommand\thesubsection {\thesection.\@arabic\c@subsection}
\renewcommand\thesubsubsection{\thesubsection.\@arabic\c@subsubsection}
\newcommand{\appendix}{
  \setcounter{section}{0}
  \renewcommand{\thesection}{\Alph{section}}}
\newcommand\abstractname{Abstract}
\newenvironment{abstract}{
  \small
  \begin{center}%
    {\bfseries \abstractname\vspace{-.5em}\vspace{\z@}}
  \end{center}\begin{quotation}}{\end{quotation}\normalsize}
