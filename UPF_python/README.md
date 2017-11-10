This is an attempt to make pseudopotential in UPF usable for ffr-LFDFT
package.
I used Python's ElementTree XML reader module.
This module works for most UPF version 2, but currently '&' sign must not be
present in the UPF file.
