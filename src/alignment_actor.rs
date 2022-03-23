use actix::prelude::*;
use std::collections::HashMap;
use crate::extractor::{align_unknown_orientation_read_u8_ref, extract_tagged_sequences, align_unknown_orientation_read};
use std::fs::File;
use actix::Recipient;
use std::io::Write;

#[derive(Message)]
#[rtype(result = "()")]
pub struct ReadRefAligner {
    pub read: String,
    pub reference: String,
    pub output_actor: Recipient<AlignedSequences>,
}

#[derive(Message)]
#[rtype(result = "()")]
pub struct AlignedSequences {
    pub read_aligned: String,
    pub ref_aligned: String,
    pub extracted_features: HashMap<String, String>,
}



pub struct Aligner;

impl Actor for Aligner {
    type Context = SyncContext<Self>;
}

impl Handler<ReadRefAligner> for Aligner {
    type Result = ();

    fn handle(&mut self, msg: ReadRefAligner, _: &mut Self::Context) -> () {
        //println!("sending stuff");
        let aligned_read1 = align_unknown_orientation_read(&msg.read, &msg.reference);
        let features = extract_tagged_sequences(&aligned_read1.2, &aligned_read1.1);


        msg.output_actor.do_send(AlignedSequences{read_aligned: aligned_read1.2,ref_aligned: aligned_read1.1,extracted_features: features});
    }
}

// our writer, which we'll have only one of
pub(crate) struct AlignmentWriter {
    pub output: File,
}

impl Actor for AlignmentWriter {
    type Context = Context<Self>;

}

// now we need to implement `Handler` on `Calculator` for the `Sum` message.
impl Handler<AlignedSequences> for AlignmentWriter {
    type Result = ();

    fn handle(&mut self, msg: AlignedSequences, ctx: &mut Context<Self>) -> () {
        //println!("Writing stuff");
        write!(self.output, ">read1\n{}\n", msg.read_aligned).unwrap();
        write!(self.output, ">ref\n{}\n", msg.ref_aligned).unwrap();
        write!(self.output, ">special\n{}\n", msg.extracted_features.iter().map(|(s,t)| format!("{}{}",&**s,&**t)).collect::<Vec<_>>().join(", ")).unwrap();
    }
}