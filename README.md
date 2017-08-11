# mmej
Microhomology predictor for microhomology-mediated end joining (MMEJ)

## Usage
```
mmej.py [-h] [--sequence-file SEQUENCE] [--good-guides GOOD_GUIDES]
             [--bad-guides BAD_GUIDES] [--max-flanking MAX_FLANKING_LENGTH]
             [--min-flanking MIN_FLANKING_LENGTH]
             [--length-weight LENGTH_WEIGHT] [--rank RANK]

Optional arguments:
  -h, --help            show this help message and exit
  --sequence-file SEQUENCE, -i SEQUENCE
                        File with the target sequence.
  --good-guides GOOD_GUIDES, -g GOOD_GUIDES
                        List of good guide sequences.
  --bad-guides BAD_GUIDES, -b BAD_GUIDES
                        List of bad guide sequences.
  --max-flanking MAX_FLANKING_LENGTH, -M MAX_FLANKING_LENGTH
                        Max length of flanking region.
  --min-flanking MIN_FLANKING_LENGTH, -m MIN_FLANKING_LENGTH
                        Min length of flanking region.
  --length-weight LENGTH_WEIGHT, -w LENGTH_WEIGHT
                        Length weight - used in scoring.
  --rank RANK           Output rangking of guides.
```
