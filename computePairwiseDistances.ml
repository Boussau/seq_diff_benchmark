(*
Installing core:
  opam install core
Compiling into bytecode:
  ocamlbuild -use-ocamlfind -pkgs core -tag thread computePairwiseDistances.byte
Compiling into native code:
  ocamlbuild -use-ocamlfind -pkgs core -tag thread computePairwiseDistances.native
*)


open Core.Std

type fasta_line =
  | Description of string
  | Partial_sequence of string

type fasta_item = {
  description : string ;
  sequence : string ;
}

type fasta_parser_state =
  | Init
  | Reading of string * string list

type alignment = fasta_item array

let fasta_line_parser l =
  if l <> "" && l.[0] = '>' then
    Description (String.slice l 1 0)
  else
    Partial_sequence l

let rev_concat xs = String.concat (List.rev xs)

let fasta_parser fn =
  In_channel.with_file fn ~f:(fun ic ->
      let rec loop i state accu =
        let input =
          In_channel.input_line ic
          |> Option.map ~f:fasta_line_parser
        in
        match state, input with
        | Init, None -> accu

        | Init, Some (Description d) ->
          loop (i + 1) (Reading (d, [])) accu

        | Init, Some (Partial_sequence _) ->
          failwithf "Malformed FASTA: expected description at line %d" i ()

        | Reading (d, []), None ->
          failwithf "Malformed FASTA: no sequence for last item %s" d ()

        | Reading (_, []), Some (Description _) ->
          failwithf "Malformed FASTA: expected sequence at line %d" i ()

        | Reading (description, seqs), Some (Description d) ->
          let sequence = rev_concat seqs in
          loop (i + 1) (Reading (d, [])) ({ description ; sequence } :: accu)

        | Reading (description, seqs), None ->
          let sequence = rev_concat seqs in
          { description ; sequence } :: accu

        | Reading (d, seqs), Some (Partial_sequence s) ->
          loop (i + 1) (Reading (d, s :: seqs)) accu
      in
      loop 0 Init []
      |> List.rev
    )


let unparse_fasta_item oc { description ; sequence } =
  fprintf oc ">%s\n" description ;
  fprintf oc "%s" sequence

let filter_map_item it =
  match String.lsplit2 it.description ~on:':' with
  | None -> None
  | Some (_, "0") -> None
  | Some (left, _) ->
    Some { it with description = left }


let seq_diffs x y =
  if String.length x.sequence <> String.length y.sequence
  then failwithf "Error: sequences %s and %s do not have the same length." x.description y.description () ;
  String.foldi x.sequence ~init:0 ~f:(fun i accu xc ->
      let yc = y.sequence.[i] in
      if xc <> yc && xc <> 'N' && yc <> 'N' then accu + 1 else accu
    )

(*
let seq_diffs xs ys =
  let n = String.length xs in
  let sum = ref 0 in
  for k = 0 to n - 1 do
    let xc = xs.[k]
    and yc = ys.[k] in
    if xc <> yc && xc <> 'N' && yc <> 'N' then incr sum
  done ;
  !sum
*)

let output_diffs oc al =
  fprintf oc "from\tto\tdist\n" ;
  let n = Array.length al in
  for i = 0 to n - 2 do
    for j = i + 1 to n - 1 do
      (*let d = seq_diffs al.(i).sequence al.(j).sequence in*)
      let d = seq_diffs al.(i) al.(j) in
      fprintf oc "%s\t%s\t%d\n" al.(i).description al.(j).description d
    done
  done


let main fa output () =
  Out_channel.with_file output ~f:(fun oc ->
      fasta_parser fa
      |> List.filter_map ~f:filter_map_item
      |> Array.of_list
      |> (fun x ->
          printf "Number of non-0 haplotypes: %d\n" (Array.length x) ;
          x)
      |> output_diffs oc
    )

let spec =
  let open Command.Spec in
  empty
  +> flag "--fasta" (required file) ~doc:"FASTA Input FASTA file"
  +> flag "--output" (required file) ~doc:"PATH Output file"

let command =
  Command.basic
    ~summary:"Computes pairwise edit distances for all pairs of sequences"
    spec
    main

let () = Command.run command
